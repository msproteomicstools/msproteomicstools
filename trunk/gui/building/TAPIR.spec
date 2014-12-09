# -*- mode: python -*-

import guidata
import guiqwt
import pymzml
import pyopenms

def dir_files(path, rel):
    ret = []
    for p,d,f in os.walk(path):
        relpath = p.replace(path, '')[1:]
        for fname in f:
            ret.append((os.path.join(rel, relpath, fname),
                        os.path.join(p, fname), 'DATA'))
    return ret

block_cipher = None


a = Analysis(['/Users/georger/ETHZ/Projects/msproteomicstools/msproteomicstools/gui/TAPIR.py'],
             pathex=['/Users/georger/Desktop'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None,
             cipher=block_cipher)
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
          name='TAPIR',
          debug=False,
          strip=None,
          upx=True,
          console=False , icon='/Users/georger/ETHZ/Projects/msproteomicstools/msproteomicstools/gui/building/icons/tapir.icns')
app = BUNDLE(exe,
             name='TAPIR.app',
             icon='/Users/georger/ETHZ/Projects/msproteomicstools/msproteomicstools/gui/building/icons/tapir.icns',
             bundle_identifier='ch.ethz.imsb.tapir',
             version='1.0',)
