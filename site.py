from app import create_app

import os
from flask import Flask

template_dir = os.path.join(os.path.dirname(__file__), 'templates')
static_dir = os.path.join(os.path.dirname(__file__), 'static')

app = Flask(__name__, template_folder=template_dir, static_folder=static_dir)

create_app(app)

if __name__ == "__main__":
    app.run(debug=True)
