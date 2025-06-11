from app import create_app

import os
from flask import Flask

base_dir = os.path.dirname(__file__)
template_dir = os.path.join(base_dir, 'templates')
static_dir = os.path.join(base_dir, 'app', 'static')

app = Flask(__name__, template_folder=template_dir, static_folder=static_dir)

create_app(app)

if __name__ == "__main__":
    app.run(debug=True)