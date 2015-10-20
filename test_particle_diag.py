from opmd_viewer.addons.pic.lpa_diagnostics import LpaDiagnostics
import pylab as plt
i = 0
select = {'z': [0, 135]}
# select = None
print('Loading Data...')
data = LpaDiagnostics(
   "/afs/desy.de/group/m/uhh/aggruener/sjalas/Simulations/opmdTest/diags/hdf5")
print('Data loaded')


mg, stdg = data.get_rms_gamma(iteration=i, species='electrons',
                              select=select)
print('RMS gamma and standard deviation: %s %s ' % (mg, stdg))

charge = data.get_charge(iteration=i, species='electrons', select=select)
print('Charge: %s C' % charge)

divx, divy = data.get_divergence(iteration=i, species='electrons',
                                 select=select)
print('Divergence: X: %s rad, Y: %s rad' % (divx, divy))

emix, emiy = data.get_emittance(iteration=i, species='electrons',
                                select=select)
print('Emittance: X: %s m rad, Y: %s m rad' % (emix, emiy))

current, zbins = data.get_current(iteration=i, species='electrons',
                                  select=select)
plt.figure()
plt.plot(zbins, current)
plt.xlabel(r'$z (\mu m)$', fontsize=18)
plt.ylabel(r'$I (A)$', fontsize=18)
plt.suptitle(r'$\mathrm{Beam \, Current}$', fontsize=22)

freq = data.get_mean_frequency(iteration=i, pol=0)
print('Mean angular Frequency: %s rad/s' % freq)

ctau = data.get_ctau(iteration=i, pol=0)
print ('Laser ctau: %s m' % ctau)

a0 = data.get_a0(iteration=i, pol=0)
print('a0: %s' % a0)

waist = data.get_laser_waist(iteration=i, pol=0)
print('Laser waist: %s' % waist)

envelope = data.get_laser_envelope(iteration=i, pol=0)
field = data.get_field(iteration=i, field='E', coord='r', theta=0)
field_slice = field[0][int( field[0].shape[0] / 2), :]
plt.figure()
plt.plot(envelope[0])
plt.plot(field_slice)
plt.ylabel(r'$E (V/m)$', fontsize=18)
plt.xlabel(r'$z (\mu m)$', fontsize=18)
plt.suptitle(r'$\mathrm{Laser \, Pulse}$', fontsize=22)
plt.figure()

env = data.get_laser_envelope(iteration=i, pol=0, index='all')
plt.imshow(env[0], aspect='auto', cmap='gist_heat',
           extent=env[1].imshow_extent*1.e6)
plt.ylabel(r'$r (\mu m)$', fontsize=18)
plt.xlabel(r'$z (\mu m)$', fontsize=18)
plt.suptitle(r'$\mathrm{Laser \, Envelope}$', fontsize=22)


wigner = data.wigner_transform(iteration=i, pol=0)
plt.figure()
plt.imshow(wigner[0], aspect='auto', cmap='bone_r', extent=wigner[1])
plt.show()
