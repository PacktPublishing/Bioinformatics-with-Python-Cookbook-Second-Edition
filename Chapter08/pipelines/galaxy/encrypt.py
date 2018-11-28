"Encrypt an YAML file with the script configuration"

import base64
import getpass
from io import StringIO
import os

from ruamel.yaml import YAML

from cryptography.fernet import Fernet
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC

password = getpass.getpass('Please enter the password:').encode()

salt = os.urandom(16)
kdf = PBKDF2HMAC(algorithm=hashes.SHA256(), length=32, salt=salt,
                 iterations=100000, backend=default_backend())
key = base64.urlsafe_b64encode(kdf.derive(password))
fernet = Fernet(key)

with open('salt', 'wb') as w:
    w.write(salt)


yaml = YAML()

content = yaml.load(open('galaxy.yaml', 'rt', encoding='utf-8'))
print(type(content), content)
output = StringIO()
yaml.dump(content, output)
print ('Encrypting:\n%s' % output.getvalue())

enc_output = fernet.encrypt(output.getvalue().encode())

with open('galaxy.yaml.enc', 'wb') as w:
    w.write(enc_output)


print("Complete, the clear version should be deleted now")
