# File: BDTReweight_api.py
# Brief: The python script to load pickled BDT reweighters that reweight
# MINERvA medium energy numu-Carbon CCQE-like events from source MC
# GENIE v2.12.6 to target MC GENIE v3.0.4 AR23.
# See https://github.com/zihaolin2000/BDTReweight
# and https://arxiv.org/abs/2510.07463.

print('Initializing python CCQELikeBDTReweighter...')

import pickle
import sys
import gc

# Append BDTReweight package path
sys.path.append('/exp/minerva/app/users/zihaolin/REWEIGHTworkdir/')

# Disable garbage collection for better performance in cpp
gc.disable()

# Load pickled reweighters from path 
# (query via category name: 0p0n, 0pNn, 1p0n, 1pNn, 2p0n, 2pNn, others)
def load_reweighter(category):
    with open(f'/exp/minerva/data/users/zihaolin/BDTReweighters/saved_reweighters_pickle/reweighter_MINERvA_ME_numuCarbon_CCQELike_GENIEv2_to_v3AR23_1mu{category}.pkl', 'rb') as f:
        reweighter = pickle.load(f)
    return reweighter

# Load all reweighters in runtime
_rw_0p0n = load_reweighter('0p0n')
_rw_0pNn = load_reweighter('0pNn')
_rw_1p0n = load_reweighter('1p0n')
_rw_1pNn = load_reweighter('1pNn')
_rw_2p0n = load_reweighter('2p0n')
_rw_2pNn = load_reweighter('2pNn')
_rw_others = load_reweighter('others')

# ================== PREDICT FUNCTIONS ==================

# reweight_variables = [
#     'total_proton_px','total_proton_py','total_proton_pz',
#     'total_proton_KE','leading_muon_py','leading_muon_pz'
# ]
def predict_weight_0p0n(features):
    w = _rw_0p0n.predict_weight_single_event(features)
    return float(w)

# reweight_variables = [
#     'leading_neutron_px', 'leading_neutron_py', 'leading_neutron_pz',
#     'total_proton_px','total_proton_py','total_proton_pz',
#     'total_proton_KE','leading_muon_py','leading_muon_pz'
# ]
def predict_weight_0pNn(features):
    w = _rw_0pNn.predict_weight_single_event(features)
    return float(w)

# reweight_variables = [
#     'leading_proton_px','leading_proton_py','leading_proton_pz',
#     'total_proton_KE','leading_muon_py','leading_muon_pz'
# ]
def predict_weight_1p0n(features):
    w = _rw_1p0n.predict_weight_single_event(features)
    return float(w)

# reweight_variables=[
#     'leading_proton_px','leading_proton_py','leading_proton_pz',
#     'total_proton_KE','leading_muon_py','leading_muon_pz',
#     'leading_neutron_px', 'leading_neutron_py', 'leading_neutron_pz'
# ]
def predict_weight_1pNn(features):
    w = _rw_1pNn.predict_weight_single_event(features)
    return float(w)

# reweight_variables = [
#     'leading_proton_px','leading_proton_py','leading_proton_pz',
#     'total_proton_KE','leading_muon_py','leading_muon_pz',
#     'subleading_proton_px', 'subleading_proton_py', 'subleading_proton_pz'
# ]
def predict_weight_2p0n(features):
    w = _rw_2p0n.predict_weight_single_event(features)
    return float(w)

# reweight_variables = [
#     'leading_proton_px','leading_proton_py','leading_proton_pz',
#     'total_proton_KE','leading_neutron_px', 'leading_neutron_py',
#     'leading_neutron_pz','leading_muon_py','leading_muon_pz',
#     'subleading_proton_px', 'subleading_proton_py', 'subleading_proton_pz'
# ]
def predict_weight_2pNn(features):
    w = _rw_2pNn.predict_weight_single_event(features)
    return float(w)

# reweight_variables = [
#     'leading_proton_px','leading_proton_py','leading_proton_pz',
#     'total_proton_KE','leading_muon_py','leading_muon_pz'
# ]
def predict_weight_others(features):
    w = _rw_others.predict_weight_single_event(features)
    return float(w)


print('Reweighter initialized.')
