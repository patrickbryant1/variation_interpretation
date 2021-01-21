
import pdb
import requests, sys

server = "http://rest.ensembl.org"
ext = "/variant_recoder/homo_sapiens"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
r = requests.post(server+ext, headers=headers, data='{ "ids" : ["A0PJY2:p.H278Y" ] }')


if not r.ok:
  r.raise_for_status()
  sys.exit()
pdb.set_trace()

decoded = r.json()
print(repr(decoded))
pdb.set_trace()
