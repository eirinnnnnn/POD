import os
ID="n128_k100_highsnr"
OP="1"*(7*64); Q="?"*128
T="""[AWGN]
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

use_AED                        = {use}
automorphism_src               = {aut}
aed_L                          = {M}

[Monte_Carlo]
iter_min                       = 0
iter_max                       = {itmax}
error_min                      = 0
error_max                      = {errmax}
monitor_slot_size              = 200000
seed_string                    = 111511015
"""
runs = {
 "SC":        dict(L=1,  aut="", M=0),
 "SCL32":     dict(L=32, aut="", M=0),
 "AE32SC_LTA":dict(L=1,  aut="./aut_lta_M32.txt",    M=32),
 "AE32SC_UTL":dict(L=1,  aut="./aut_utl_M32.txt",    M=32),
 "AE32SC_PU": dict(L=1,  aut="./aut_pu_M32.txt",     M=32),
 "AE32SC_RND":dict(L=1,  aut="./aut_random_M32.txt", M=32),
 "AE4SCL8_UTL":dict(L=8, aut="./aut_utl_M4.txt",     M=4),
}
import sys
start,end,itmax,errmax,suffix = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],(sys.argv[5] if len(sys.argv)>5 else "")
for name,cfg in runs.items():
    open(name+suffix+".ini","w").write(T.format(
        start=start,end=end,op=OP,q=Q,id=ID,L=cfg["L"],
        use=("true" if cfg["aut"] else "false"),aut=cfg["aut"],M=cfg["M"],
        itmax=itmax,errmax=errmax))
print("wrote", len(runs), "inis")
