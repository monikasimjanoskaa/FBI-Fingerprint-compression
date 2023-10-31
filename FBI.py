import numpy as np
import matplotlib 
from matplotlib import pyplot as plt
import pywt
#%% Vchituvanje na slikata
img=matplotlib.image.imread('fingerprint.tif')

#pretvaranje na vrednostite na matricata od 16-bitni broevi vo decimalni  broevi vo opseg od 0 do 1
def im2double(im):
    info=np.iinfo(im.dtype)               
    return im.astype(np.float) / info.max
f0=im2double(img)

#prikazuvanje na slikata vo format grayscale
plt.imshow(f0, cmap='gray')
plt.title('Originalna slika')

#%% Definiranje na mek i tvrd prag
def soft_treshold(signal,level):
 temp=(abs(signal)-level) 
 temp=(temp+abs(temp))/2 
 y=np.multiply(np.sign(signal),temp)
 return y

def hard_treshold(signal,level):
 y=np.multiply(signal,abs(signal)>level) 
 return y
#%% LEVEL 2
#paketska vejvlet transformacija vrz signalot
coeff=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=2)

#bidekji koeficientite gi dobivame vo lista,potrebno e da gi izdvoime vo soodvetni nizi
jazol=[node.path for node in coeff.get_level(2,'natural')]

#presmetka na vkupniot broj na nenulti koeficienti vo listata
# bidejki imame izbrano level 2 dekompozicija dobivame vkupno 16 nizi
vkbroj = 0
for i in range (0,16):
    vkbroj = vkbroj + np.count_nonzero(coeff[jazol[i]].data)

#%% Soft treshold so prag = 0.01

temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=2) #vo temp gi zacuvuvame originalnite koeficienti
jazol=[node.path for node in temp.get_level(2,'natural')] 

#dejstvuvanje na pragot za koeficientite na detali
prag=0.01
for i in range(1,16): 
    temp[jazol[i]].data=soft_treshold(temp[jazol[i]].data,prag)
    
#rekonstrukcija na slikata
f2=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f2, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level 2 ,prag=0.01 ,soft')

#presmetuvanje na brojot na neluti elementi posle dejstvuvanje na pragot
vkbroj_temp=0
for i in range (0,16):
    vkbroj_temp = vkbroj_temp + np.count_nonzero(temp[jazol[i]].data)

#PRD
#f0-originalna slika 
#f2-rekonstruirana slika 
x1=(f0-f2)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

#PSNR
#f0-originalna slika 
#f2-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f2,'fro'))**2/(np.size(f0))
PSNR=10*np.log10((f0.max())**2/MSE)

#CR                                                 
CR=vkbroj/vkbroj_temp

#%% Hard treshold so prag = 0.01

temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=2) #vo temp gi zacuvuvame originalnite koeficienti
jazol=[node.path for node in temp.get_level(2,'natural')] 

#dejstvuvanje na pragot samo za koeficientite na detali
prag = 0.01
for i in range(1,16):
    temp[jazol[i]].data=hard_treshold(temp[jazol[i]].data,prag)
    
#rekonstrukcija
f2=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f2, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level2,prag=0.01,hard')
    
#presmetuvanje na brojot na neluti elementi posle dejstvuvanje na pragot
vkbroj_temp=0
for i in range (0,16):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol[i]].data)

#PRD
#f0-originalna slika 
#f2-rekonstruirana slika 
x1=(f0-f2)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

#PSNR
#f0-originalna slika 
#f2-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f2,'fro'))**2/(np.size(f0))
PSNR=10*np.log10((f0.max())**2/MSE)

#CR
CR=vkbroj/vkbroj_temp

#%% Soft treshold so prag=0.7
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=2) #vo temp gi zacuvuvame originalnite koeficienti
jazol=[node.path for node in temp.get_level(2,'natural')] 

#dejstvuvanje na pragot samo za koeficientite na detali
prag=0.7
for i in range(1,16): 
    temp[jazol[i]].data=soft_treshold(temp[jazol[i]].data,prag)
    
#rekonstrukcija
f2=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f2, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level2 ,prag=0.7,soft)')

#presmetuvanje na brojot na neluti elementi posle dejstvuvanje na pragot
vkbroj_temp=0
for i in range (0,16):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol[i]].data)

#PRD
#f0-originalna slika 
#f2-rekonstruirana slika 
x1=(f0-f2)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

#PSNR
#f0-originalna slika 
#f2-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f2,'fro'))**2/(np.size(f0))
PSNR=10*np.log10((f0.max())**2/MSE)

#CR   
CR=vkbroj/vkbroj_temp

#%% Hard treshold so prag=0.7
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=2) #vo temp gi zacuvuvame originalnite koeficienti
jazol=[node.path for node in temp.get_level(2,'natural')] 

#dejstvuvanje na pragot samo za koeficientite na detali
prag = 0.7
for i in range(1,16):
    temp[jazol[i]].data=hard_treshold(temp[jazol[i]].data,prag)
 
#rekonstrukcija
f2=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f2, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level2,prag=0.7,hard')
    
#presmetuvanje na brojot na neluti elementi posle dejstvuvanje na pragot
vkbroj_temp=0
for i in range (0,16):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol[i]].data)
  
#PRD
#f0-originalna slika 
#f2-rekonstruirana slika 
x1=(f0-f2)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

# PSNR
#f0-originalna slika 
#f2-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f2,'fro'))**2/(np.size(f0)) 
PSNR=10*np.log10((f0.max())**2/MSE)

#CR
CR=vkbroj/vkbroj_temp

#%% Level 4
coeff=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=4)

jazol4=[node.path for node in coeff.get_level(4,'natural')]

vkbroj=0
for i in range (0,256):
    vkbroj=vkbroj + np.count_nonzero(coeff[jazol4[i]].data)

#%% Soft treshold so prag = 0.01
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=4)
jazol4=[node.path for node in temp.get_level(4,'natural')] 

prag=0.01
for i in range(1,256): 
    temp[jazol4[i]].data=soft_treshold(temp[jazol4[i]].data,prag)

#rekonstrukcija
f4=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f4, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level 4,prag=0.01,soft')

vkbroj=0
for i in range (0,256):
    vkbroj=vkbroj + np.count_nonzero(coeff[jazol4[i]].data)

vkbroj_temp=0
for i in range (0,256):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol4[i]].data)

#PRD
#f0-originalna slika 
#f4-rekonstruirana slika 
x1=(f0-f4)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

# PSNR
#f0-originalna slika 
#f4-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f4,'fro'))**2/(np.size(f0))
PSNR=10*np.log10((f0.max())**2/MSE)

#CR    
CR=vkbroj/vkbroj_temp
#%% Hard treshold so prag = 0.01
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=4)
jazol4=[node.path for node in temp.get_level(4,'natural')] 

for i in range(1,256):
    temp[jazol4[i]].data=hard_treshold(temp[jazol4[i]].data,prag)
    
#rekonstrukcija
f4=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f4, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level 4,prag=0.01,hard')

vkbroj_temp=0
for i in range (0,256):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol4[i]].data)

#PRD
#f0-originalna slika 
#f4-rekonstruirana slika 
x1=(f0-f4)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

# PSNR
#f0-originalna slika 
#f4-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f4,'fro'))**2/(np.size(f0))
PSNR=10*np.log10((f0.max())**2/MSE)


#CR
CR=vkbroj/vkbroj_temp

#%%  Soft treshold so prag = 0.7
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=4)
jazol4=[node.path for node in temp.get_level(4,'natural')] 

prag=0.7 
for i in range(1,256): 
    temp[jazol4[i]].data=soft_treshold(temp[jazol4[i]].data,prag)
    
#rekonstrukcija
f4=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f4, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level 4 ,prag=0.7,soft')

vkbroj_temp=0
for i in range (0,256):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol4[i]].data)
 
#PRD
#f0-originalna slika 
#f4-rekonstruirana slika 
x1=(f0-f4)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

# PSNR
#f0-originalna slika 
#f4-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f4,'fro'))**2/(np.size(f0))
PSNR=10*np.log10((f0.max())**2/MSE)

#CR
CR=vkbroj/vkbroj_temp

#%% Hard treshold so prag = 0.7
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=4)
jazol4=[node.path for node in temp.get_level(4,'natural')] 

for i in range(1,256):
    temp[jazol4[i]].data=hard_treshold(temp[jazol4[i]].data,prag)
    
#rekonstrukcija
f4=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f4, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level 4,prag=0.7,hard')
    
vkbroj_temp=0
for i in range (0,256):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol4[i]].data)
    
#PRD
#f0-originalna slika 
#f4-rekonstruirana slika
x1=(f0-f4)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())
  
# PSNR
#f0-originalna slika 
#f4-rekonstruirana slika
MSE=(np.linalg.norm(f0-f4,'fro'))**2/(np.size(f0))
PSNR=10*np.log10((f0.max())**2/MSE)

#CR
CR=vkbroj/vkbroj_temp


#%% LEVEL 6
coeff=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=6)

jazol=[node.path for node in coeff.get_level(6,'natural')]

# bidejki imame izbrano level 6 dekompozicija dobivame vkupno 4096 nizi 
vkbroj=0
for i in range (0,4096):
    vkbroj = vkbroj + np.count_nonzero(coeff[jazol[i]].data)

#%% Soft treshold so prag = 0.7

temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=6) #vo temp gi zacuvuvame originalnite koeficienti
jazol6=[node.path for node in temp.get_level(6,'natural')] 

#dejstvuvanje na pragot samo za koeficientite na detali
prag=0.7

for i in range(1,4096):
    temp[jazol6[i]].data=soft_treshold(temp[jazol6[i]].data,prag)
    
#rekonstrukcija
f6=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f6, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level6,prag=0.7,soft')

#presmetuvanje na brojot na neluti elementi posle dejstvuvanje na pragot
vkbroj_temp=0
for i in range (0,4096):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol6[i]].data)

#PRD
#f0-originalna slika 
#f6-rekonstruirana slika 
x1=(f0-f6)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

#PSNR
#f0-originalna slika 
#f6-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f6,'fro'))**2/(np.size(f0)) 
PSNR=10*np.log10((f0.max())**2/MSE)

#CR                                                
CR=vkbroj/vkbroj_temp

#%% Hard treshold so prag = 0.7
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=6)
jazol6=[node.path for node in temp.get_level(6,'natural')] 

prag = 0.7
for i in range(1,4096):
    temp[jazol6[i]].data=hard_treshold(temp[jazol6[i]].data,prag)
    
#rekonstrukcija
f6=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f6, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level6,prag=0.7,hard')
    
vkbroj_temp=0
for i in range (0,4096):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol6[i]].data)

#PRD
#f0-originalna slika 
#f6-rekonstruirana slika 
x1=(f0-f6)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

# PSNR
#f0-originalna slika 
#f6-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f6,'fro'))**2/(np.size(f0)) 
PSNR=10*np.log10((f0.max())**2/MSE)

# CR
CR=vkbroj/vkbroj_temp

#%% Soft treshold prag = 0.01
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=6)
jazol6=[node.path for node in temp.get_level(6,'natural')] 

prag=0.01
for i in range(1,4096):
    temp[jazol6[i]].data=soft_treshold(temp[jazol6[i]].data,prag)
    
#rekonstrukcija
f6=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f6, cmap='gray') 
plt.title('Рekonstrukcija na slikata so level 6,prag=0.01,soft')

vkbroj_temp=0
for i in range (0,4096):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol6[i]].data)
 
#PRD
#f0-originalna slika 
#f6-rekonstruirana slika 
x1=(f0-f6)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

# PSNR
#f0-originalna slika 
#f6-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f6,'fro'))**2/(np.size(f0)) 
PSNR=10*np.log10((f0.max())**2/MSE)


#CR
CR=vkbroj/vkbroj_temp

#%%  Hard treshold so prag = 0.01
temp=pywt.WaveletPacket2D(f0,'bior4.4', maxlevel=6)
jazol6=[node.path for node in temp.get_level(6,'natural')] 

prag=0.01
for i in range(1,4096):
    temp[jazol6[i]].data=hard_treshold(temp[jazol6[i]].data,prag)
 
#rekonstrukcija
f3=temp.reconstruct(update=True)
plt.figure()
plt.imshow(f6, cmap='gray') 
plt.title('Rekonstrukcija na slikata so level6, prag=0.01,hard')
    
vkbroj_temp=0
for i in range (0,4096):
    vkbroj_temp= vkbroj_temp + np.count_nonzero(temp[jazol6[i]].data)

#PRD
#f0-originalna slika 
#f6-rekonstruirana slika 
x1=(f0-f6)**2
x2=f0**2
PRD=np.sqrt(x1.sum()/x2.sum())

#PSNR
#f0-originalna slika 
#f6-rekonstruirana slika 
MSE=(np.linalg.norm(f0-f3,'fro'))**2/(np.size(f0)) 
PSNR=10*np.log10((f0.max())**2/MSE)

#CR  
CR=vkbroj/vkbroj_temp


