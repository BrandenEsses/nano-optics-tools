from flask import Flask, render_template, request, redirect, session, Blueprint
from app import app

calculators_blueprint = Blueprint('calculators', __name__, url_prefix='/calculators', template_folder="templates")

@calculators_blueprint.route('/')
def calculators():
    return render_template("calculators.html", title="Calculators")