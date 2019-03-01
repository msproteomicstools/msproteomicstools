# -*- mode: python -*-

# Important: 
# Run this from the main directory:
#   copy gui\building\TAPIR_win.spec .
#   pyinstaller.exe TAPIR_win.spec

import guidata
import guiqwt
import pymzml
import pyopenms
import msproteomicstoolslib
import PyMSNumpress
import PyQt4

import os
cwd = os.getcwd()

def dir_files(path, rel):
    ret = []
    for p,d,f in os.walk(path):
        relpath = p.replace(path, '')[1:]
        for fname in f:
            ret.append((os.path.join(rel, relpath, fname),
                        os.path.join(p, fname), 'DATA'))
    return ret

block_cipher = None

a = Analysis(['gui\\TAPIR.py'],
             pathex=[cwd],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None,
             cipher=block_cipher)

# This allows us to get data which will work when using relative paths inside
# our code (best works with:
## root = os.path.dirname(__file__)
#
a.datas.extend(dir_files(os.path.join(os.path.dirname(guidata.__file__),
    'images'), os.path.join('guidata', 'images')))
a.datas.extend(dir_files(os.path.join(os.path.dirname(guiqwt.__file__),
    'images'), os.path.join('guiqwt', 'images')))
a.datas.extend(dir_files(os.path.join(os.path.dirname(pymzml.__file__),
    'obo'), os.path.join('pymzml', 'obo')))
a.datas.extend(dir_files(os.path.join(os.path.dirname(pyopenms.__file__),
    'share'), os.path.join('pyopenms', 'share')))
pyz = PYZ(a.pure,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='TAPIR.exe',
          debug=False,
          strip=None,
          upx=True,
          console=False, # set to True for debug
          icon='gui\\building\\icons\\tapir.ico')
