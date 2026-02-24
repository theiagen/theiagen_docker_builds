FROM python:3.12-slim

RUN pip install --no-cache-dir pandas==2.3.1 requests==2.32.3

COPY freyja_to_long.py /usr/local/bin/freyja_to_long.py
COPY freyja_microreact.py /usr/local/bin/freyja_microreact.py
RUN chmod +x /usr/local/bin/freyja_to_long.py /usr/local/bin/freyja_microreact.py

ENTRYPOINT ["python", "/usr/local/bin/freyja_to_long.py"]
