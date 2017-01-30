#! /usr/bin/python3

import requests
import re
import boto
import boto.s3.connection
import shutil
import os
import tarfile
import subprocess
from boto.s3.key import Key


#set security keys
AWS_ACCESS_KEY_ID = 'AKIAJOIGNLI4C3YVBTRA'
AWS_SECRET_ACCESS_KEY = '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'
TMPDIR='/home/ubuntu/tmp'
BACHOME='/home/ubuntu/BacBuilder'


#get the bucket id
r = requests.get("http://169.254.169.254/latest/user-data")

if r.status_code != 200:
    exit(-1)

m=re.search('BUCKETID: (.*)', r.text)

#check we have a bucket id, else fail
if m:
    bucketid = m.group(1)

else:
    print ("No data bucket specified, exiting!")
    exit(-1);

#Check if we need to drop to admin mode
if re.match('ADMIN', bucketid):
    print ("Entering admin mode")
    exit(0)


print ("Checking user data")

#download the input tar and param list
print ("Downloading input from "+bucketid+" bucket")


conn = boto.connect_s3(AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY)
bucket = conn.get_bucket(bucketid)

modeltarball = bucket.get_key('model.tgz')

if modeltarball is None:
    print ("Couldn't download model from bucket "+bucketid)
    exit (-1)

shutil.rmtree(TMPDIR, ignore_errors=True)
os.mkdir(TMPDIR)

modeltarball.get_contents_to_filename(TMPDIR+'/model.tgz')

#unpack the tar
tar = tarfile.open(TMPDIR+'/model.tgz')
print("Extracting:")
tar.list()
tar.extractall(path=TMPDIR)
tar.close()

#move the files to the bac directories
files = os.listdir(TMPDIR)

pdb=""

for file in files:
    if re.match('.*.pdb', file):
        m=re.search('(.*).pdb', file)
        if m:
            pdb=m.group(1)
            print("PDB: "+pdb)
            pdbdir=BACHOME+"/raw_pdbs/egfr/"
            shutil.rmtree(pdbdir, ignore_errors=True)
            os.mkdir(pdbdir)
            print("Moving file "+file+" to "+pdbdir)
            shutil.move(TMPDIR+"/"+file, pdbdir)
        else:
            exit(-1)
    elif os.path.isdir(TMPDIR+"/"+file):
        drugdir=BACHOME+"/drugs_par/egfr/resp/"
        shutil.rmtree(drugdir, ignore_errors=True)
        os.mkdir(drugdir)
        print("Moving file "+file+" to "+drugdir)
        shutil.move(TMPDIR+"/"+file, drugdir)

#ready to run bac!!!!
print("PDB: "+pdb)

p = subprocess.run([BACHOME+'/exe/builder.pl', '-enzyme', 'egfr', '-pdb', pdb, '2>&1'], stdout=subprocess.PIPE)

print ("return code: "+str(p.returncode))

stdout=str(p.stdout, encoding='UTF-8')

print ("Stdout: "+stdout)

m=re.search('creating.system.concourse..(.*)', stdout)
if m:
    directory=m.group(1)
    print("Output written to: "+directory)

if os.path.isfile(directory+"/build/leap.log"):
    print("========LEAP LOG========")
    file = open(directory+"/build/leap.log","r")
    print(file.read())

if os.path.isfile(directory+"/build/build.log"):
    print("========BUILD LOG========")
    file = open(directory+"/build/build.log","r")
    print(file.read())

print("=======TARRING=======")
tar2 = tarfile.open(TMPDIR+'/built.tgz', "w:gz")
tar2.add(directory, arcname='rep0')
print("Compressed:")
tar2.list()
tar2.close()

buildtar = Key(bucket)

buildtar.key='build.tgz'
buildtar.set_contents_from_filename(TMPDIR+'/built.tgz')

print("Uploaded built.tgz to "+bucketid)

exit(-1)
