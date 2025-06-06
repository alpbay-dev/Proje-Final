def create_app(app):
    from .routes import main
    app.register_blueprint(main)
    return app