FROM python:3.8.1-buster
COPY Pipfile Pipfile
COPY Pipfile.lock Pipfile.lock

WORKDIR /app
RUN pip install --upgrade pip
RUN pip install pipenv
RUN pipenv install --system
EXPOSE 8080
COPY . /app
CMD python server.py

