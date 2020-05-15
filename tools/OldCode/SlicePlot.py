import yt as yt
ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

A = ['HI density','HII density','HeIII density', 'HeII density','metallicity','pressure','temperature']
for a in A:
    slc = yt.SlicePlot(ds, 'z',a)
    slc.save('Slice_'+a+'_z_rei20c1_a0.1667')
print('Done')
