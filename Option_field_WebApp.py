from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField, PasswordField
from wtforms.validators import DataRequired, Length

# https://werkzeug.palletsprojects.com/en/1.0.x/datastructures/#werkzeug.datastructures.FileStorage
# https://wtforms.readthedocs.io/en/latest/fields/#wtforms.fields.FileField
# https://stackoverflow.com/questions/20015550/read-file-data-without-saving-it-in-flask
# http://izmailoff.github.io/web/flask-file-streaming/
# https://pythonise.com/series/learning-flask/flask-uploading-files

# Text Area:
# https://stackoverflow.com/questions/7979548/how-to-render-my-textarea-with-wtforms
# https://stackoverflow.com/questions/12099741/how-do-you-set-a-default-value-for-a-wtforms-selectfield


class Login(FlaskForm):
    """ Login form"""

    # for fastq files: min and max length of reads
    # file upload
    name = StringField('Username', validators=[DataRequired(), Length(min=3, max=15)])
    password = PasswordField('Password', validators=[DataRequired()])
    submit = SubmitField('Login')


