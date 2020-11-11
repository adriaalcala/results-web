FROM python:3.8.1-buster
COPY . /app
WORKDIR /app
RUN pip install --upgrade pip
RUN pip install pipenv
RUN pipenv install
EXPOSE 8080
CMD pipenv run python server.py