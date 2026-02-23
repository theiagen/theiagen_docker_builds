FROM python:3.12-slim

RUN pip install --no-cache-dir pandas==2.3.1

COPY freyja_to_long.py /usr/local/bin/freyja_to_long.py

ENTRYPOINT ["python", "/usr/local/bin/freyja_to_long.py"]
