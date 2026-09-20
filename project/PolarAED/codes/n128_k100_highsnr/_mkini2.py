import sys
ID="n128_k100_highsnr"
OP="1"*(7*64); Q="?"*128
AE="""[AWGN]
start                          = {start}
end                            = {end}
step                           = 0.5
step_type                      = SNR
seed_string                    = -2

[AdjustPolarDecoder]
matrix_src                     = byHmatrix
Gmatrix_path                   =
Hmatrix_path                   = ./{id}_H.matrix
operationArray                 = {op}
bha_value_setting              = {q}
target_raw_BER                 = 0.01
list_size                      = {L}
permutation_src                =
permutation_random_seed        = -1
OnlyInit                       = false

use_AED                        = true
automorphism_src               = ./aut_utl_M{M}.txt
aed_L                          = {M}

[Monte_Carlo]
iter_min                       = 0
iter_max                       = {itmax}
error_min                      = 0
error_max                      = {errmax}
monitor_slot_size              = 200000
seed_string                    = 111511015
"""
OSD="""[AWGN]
start                          = {start}
end                            = {end}
step                           = 0.5
step_type                      = SNR
seed_string                    = -2

[OSD]
matrix_src                     = byGmatrix
Gmatrix_path                   = ./{id}.matrix
Hmatrix_path                   =
OSD_order                      = {order}

[Monte_Carlo]
iter_min                       = 0
iter_max                       = {itmax}
error_min                      = 0
error_max                      = {errmax}
monitor_slot_size              = 20000
seed_string                    = 111511015
"""
start,end,itmax,errmax = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
sfx = sys.argv[5] if len(sys.argv)>5 else ""
# effective list size M * L = 32 for every AE row
for name,(M,L) in {"AE16SCL2_UTL":(16,2),"AE8SCL4_UTL":(8,4),"AE2SCL16_UTL":(2,16)}.items():
    open(name+sfx+".ini","w").write(AE.format(start=start,end=end,op=OP,q=Q,id=ID,
                                              L=L,M=M,itmax=itmax,errmax=errmax))
open("OSD1"+sfx+".ini","w").write(OSD.format(start=start,end=end,id=ID,order=1,
                                             itmax=itmax,errmax=errmax))
print("wrote 4 inis")
