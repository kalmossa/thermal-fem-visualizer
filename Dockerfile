# Image de production du visualiseur FEM.
#
# Le serveur de développement Flask n'est pas fait pour être exposé : le service
# passe par gunicorn, et `debug` reste à False dans app.py — le débogueur
# Werkzeug permet l'exécution de code arbitraire à distance.
FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1     PYTHONDONTWRITEBYTECODE=1     PORT=8000

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY fem_ballast_beton.py materials.py track_beam.py app.py fem_visualizer.html ./

# Un utilisateur non privilégié : le service n'a besoin d'écrire nulle part.
RUN useradd --create-home --shell /usr/sbin/nologin fem && chown -R fem:fem /app
USER fem

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=4s --start-period=10s   CMD python -c "import urllib.request,os,sys; sys.exit(0 if urllib.request.urlopen(f'http://127.0.0.1:{os.environ.get(\"PORT\",8000)}/api/health', timeout=3).status == 200 else 1)"

# Deux workers : un calcul mobilise un cœur pendant ~0,7 s et app.py sérialise
# déjà les calculs au sein d'un worker. Le timeout couvre le maillage le plus
# fin autorisé.
CMD gunicorn --bind 0.0.0.0:$PORT --workers 2 --threads 4 --timeout 60 app:app
