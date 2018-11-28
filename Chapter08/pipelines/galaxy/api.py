import base64
from collections import defaultdict
import ftplib
import getpass
import pprint
import warnings

from ruamel.yaml import YAML

from cryptography.fernet import Fernet
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC

import pandas as pd

from bioblend.galaxy import GalaxyInstance

pp = pprint.PrettyPrinter()
warnings.filterwarnings('ignore')
# explain above, and warn


with open('galaxy.yaml.enc', 'rb') as f:
    enc_conf = f.read()


password = getpass.getpass('Please enter the password:').encode()
with open('salt', 'rb') as f:
    salt = f.read()
kdf = PBKDF2HMAC(algorithm=hashes.SHA256(), length=32, salt=salt,
                 iterations=100000, backend=default_backend())
key = base64.urlsafe_b64encode(kdf.derive(password))
fernet = Fernet(key)

yaml = YAML()
conf = yaml.load(fernet.decrypt(enc_conf).decode())

server = conf['server']
rest_protocol = conf['rest_protocol']
rest_port = conf['rest_port']
user = conf['user']
password = conf['password']
ftp_port = int(conf['ftp_port'])
api_key = conf['api_key']

rest_url = '%s://%s:%d' % (rest_protocol, server, rest_port)

history_name = 'bioinf_example'

gi = GalaxyInstance(url=rest_url, key=api_key)
gi.verify = False
histories = gi.histories

print('Existing histories:')
for history in histories.get_histories():
    if history['name'] == history_name:
        histories.delete_history(history['id'])
    print('  - ' + history['name'])
print()

ds_history = histories.create_history(history_name)


print('Uploading file')
ftp = ftplib.FTP() 
ftp.connect(host=server,port=ftp_port)
ftp.login(user=user, passwd=password)
f = open('LCT.bed','rb')
#ftp.set_pasv(False)  # explain
ftp.storbinary('STOR LCT.bed', f)
f.close() 
ftp.close()

gi.tools.upload_from_ftp('LCT.bed', ds_history['id'])
print()

contents = gi.histories.show_history(ds_history['id'], contents=True)

def summarize_contents(contents):
    summary = defaultdict(list)
    for item in contents:
        summary['íd'].append(item['id'])
        summary['híd'].append(item['hid'])
        summary['name'].append(item['name'])
        summary['type'].append(item['type'])
        summary['extension'].append(item['extension'])
    return pd.DataFrame.from_dict(summary)

print('History contents:')
pd_contents = summarize_contents(contents)
print(pd_contents)
print()

print('Metadata for LCT.bed')
bed_ds = contents[0]
pp.pprint(bed_ds)
print()

print('Metadata about all tools')
all_tools = gi.tools.get_tools()
pp.pprint(all_tools)
print()

bed2gff = gi.tools.get_tools(name='Convert BED to GFF')[0]
print("Convert BED to GFF metadata:")
pp.pprint(gi.tools.show_tool(bed2gff['id'], io_details=True, link_details=True))
print()

def dataset_to_param(dataset):
    return dict(src='hda', id=dataset['id'])

tool_inputs = {
    'input1': dataset_to_param(bed_ds)
    }

#hid!


gi.tools.run_tool(ds_history['id'], bed2gff['id'], tool_inputs=tool_inputs)
