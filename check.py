import os

count = 0
dft_count = 0
for mp_id  in os.listdir('/home/jgrandel/mace_phonons/dft_pbe_phonons/'):
    count += 1
    if os.path.exists(f'/home/jgrandel/mace_phonons/phonons_mp_0b3_medium/{mp_id}'):
        dft_count += 1
    else:
        print(mp_id)
print(count)
print(dft_count)