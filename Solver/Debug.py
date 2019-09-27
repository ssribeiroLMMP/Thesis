import matplotlib.pyplot as plt
import matplotlib.cm as cm

mycmap = cm.get_cmap('jet')
impath = '/home/sergio/Documents/HeleShawCode/FCS_v1.0/PostProcessing/Cases/vargesPR_HeleShawCell_SurfTensTest_sigma0_01_miStar_2e-1/Images/'

DD = VectorFunctionSpace(meshObj,'CG',1)
DG = VectorFunctionSpace(meshObj,'CG',2)

plt.figure(num=1, figsize=(30, 15), dpi=180, facecolor='w', edgecolor='k')
plt.clf
cax = plot(c0)
plt.colorbar(cax,orientation='horizontal', cmap = mycmap)
plt.savefig(impath+'C0Test.png')

plt.figure(num=2, figsize=(30, 15), dpi=180, facecolor='w', edgecolor='k')
plt.clf
cax = plot(c1)
plt.colorbar(cax,orientation='horizontal', cmap = mycmap)
plt.savefig(impath+'C1Test.png')

plt.figure(num=3, figsize=(30, 15), dpi=180, facecolor='w', edgecolor='k')
plt.clf
cax = plot(c0)
plt.colorbar(cax,orientation='horizontal', cmap = mycmap)
plt.savefig(impath+'CReinitTest.png')

plt.figure(num=2, figsize=(30, 15), dpi=180, facecolor='w', edgecolor='k')
plot(c0)
plt.savefig(impath+'CReinitTest.png')

plt.figure(num=1, figsize=(30, 15), dpi=180, facecolor='w', edgecolor='k')

fId = open(impath+'C0Values.txt','w+')
# Coordinates
x = meshObj.coordinates()[:,0]
y = meshObj.coordinates()[:,1]

for i in range(0,len(c0.vector())):
    fId.write('x = '+str(x[i])+' , y = '+str(y[i])+' : c = ')
    fId.write(str(c0.vector()[i]))
    fId.write('\n')

fId.close