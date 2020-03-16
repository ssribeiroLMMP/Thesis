# Debug Routines
import dolfin as df
import matplotlib.pyplot as plt

# # Open New Figure
# fig,axs=plt.subplots(2,2)
# # Plot Interface Variables
# plt.axes(axs[0,0]); df.plot(c0); axs[0,0].set_title('Concentration Field')
# plt.axes(axs[1,0]); df.plot(grad_c); axs[1,0].set_title('Gradient of Concentration')
# plt.axes(axs[0,1]); df.plot(nGamma); axs[0,1].set_title('Normal Vectors')
# plt.axes(axs[1,1]); df.plot(deltaDirac); axs[1,1].set_title('DeltaDirac')
# plt.show()


# # Open New Figure
fig,axs=plt.subplots(2,3)
# Plot Interface Tension Term
plt.axes(axs[0,0]); df.plot(rho); axs[0,0].set_title('Density Field')
plt.axes(axs[1,0]); df.plot(grad_rhoNorm); axs[1,0].set_title('Normalized Density Gradient')
plt.axes(axs[0,1]); df.plot(rhoNorm); axs[0,1].set_title('Normalized Density')
plt.axes(axs[1,1]); df.plot(dDirac); axs[1,1].set_title('Dirac Delta')
plt.axes(axs[0,2]); df.plot(k); axs[0,2].set_title('Curvature')
plt.axes(axs[1,2]); df.plot(teste); axs[1,2].set_title('Surface Tension Term')
plt.show()

# # Open New Figure
# fig,axs=plt.subplots()
# # Plot Interface Tension Term
# plt.axes(axs); df.plot(c0); axs.set_title('Concentrations')
# plt.show()