from waitress import serve
from app import server


if __name__ == "__main__":
    serve(server)