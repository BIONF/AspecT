import pickle
from flask_bcrypt import Bcrypt

bcrypt = Bcrypt()

user = input("Please Enter the new Username: ")
while True:
    user2 = input("Please confirm the new Username: ")
    if user == user2:
        break
    else:
        print('Please make sure that the usernames match!')

pw = input("Please Enter the new password: ")
while True:
    pw2 = input("Please confirm the new password: ")
    if pw == pw2:
        break
    else:
        print('Please make sure that the passwords match!')

ls = [bcrypt.generate_password_hash(user), bcrypt.generate_password_hash(pw)]
with open(r'config/login.txt', 'wb') as fp:
    pickle.dump(ls, fp)
