#!/bin/env python
from pomegranate import *
from Bio import AlignIO
import numpy as np

# hap1: { hap2: { ref: ...
observations = {
    "3": { 
        "A" : { "A": { "A": 0,"T": 2,"C": 2,"G": 2,"N": 14,"-": 10},
                "T": { "A": 3,"T": 1,"C": 4,"G": 4,"N": 14,"-": 11},
                "C": { "A": 3,"T": 4,"C": 1,"G": 4,"N": 14,"-": 11},
                "G": { "A": 3,"T": 4,"C": 4,"G": 1,"N": 14,"-": 11},
                "N": { "A": 14,"T": 14,"C": 14,"G": 14,"N": 14,"-": 14},
                "-": { "A": 7,"T": 8,"C": 8,"G": 8,"N": 14,"-": 5}
                },
        "T" : { "T": { "T": 0,"C": 2,"G": 2,"A": 2,"N": 14,"-": 10},
                "C": { "T": 3,"C": 1,"G": 4,"A": 4,"N": 14,"-": 11},
                "G": { "T": 3,"C": 4,"G": 1,"A": 4,"N": 14,"-": 11},
                "A": { "T": 3,"C": 4,"G": 4,"A": 1,"N": 14,"-": 11},
                "N": { "T": 14,"C": 14,"G": 14,"A": 14,"N": 14,"-": 14},
                "-": { "T": 7,"C": 8,"G": 8,"A": 8,"N": 14,"-": 5}
                },
        "C" : { "C": { "C": 0,"G": 2,"A": 2,"T": 2,"N": 14,"-": 10},
                "G": { "C": 3,"G": 1,"A": 4,"T": 4,"N": 14,"-": 11},
                "A": { "C": 3,"G": 4,"A": 1,"T": 4,"N": 14,"-": 11},
                "T": { "C": 3,"G": 4,"A": 4,"T": 1,"N": 14,"-": 11},
                "N": { "C": 14,"G": 14,"A": 14,"T": 14,"N": 14,"-": 14},
                "-": { "C": 7,"G": 8,"A": 8,"T": 8,"N": 14,"-": 5}
                },
        "G" : { "G": { "G": 0,"A": 2,"T": 2,"C": 2,"N": 14,"-": 10},
                "A": { "G": 3,"A": 1,"T": 4,"C": 4,"N": 14,"-": 11},
                "T": { "G": 3,"A": 4,"T": 1,"C": 4,"N": 14,"-": 11},
                "C": { "G": 3,"A": 4,"T": 4,"C": 1,"N": 14,"-": 11},
                "N": { "G": 14,"A": 14,"T": 14,"C": 14,"N": 14,"-": 14},
                "-": { "G": 7,"A": 8,"T": 8,"C": 8,"N": 14,"-": 5}
                },
        "N" : { "G": { "G": 14,"A": 14,"T": 14,"C": 14,"N": 14,"-": 14},
                "A": { "G": 14,"A": 14,"T": 14,"C": 14,"N": 14,"-": 14},
                "T": { "G": 14,"A": 14,"T": 14,"C": 14,"N": 14,"-": 14},
                "C": { "G": 14,"A": 14,"T": 14,"C": 14,"N": 14,"-": 14},
                "N": { "G": 14,"A": 14,"T": 14,"C": 14,"N": 14,"-": 14},
                "-": { "G": 14,"A": 14,"T": 14,"C": 14,"N": 14,"-": 14}
                },
        "-" : { "-": { "-": None,"A": 12,"T": 12,"C": 12,"G": 12,"N": 14},
                "A": { "-": 9,"A": 6, "T": 13, "C": 13,"G": 13,"N": 14},
                "T": { "-": 9,"A": 13,"T": 6,"C": 13,"G": 13,"N": 14},
                "C": { "-": 9,"A": 13,"T": 13,"C": 6,"G": 13,"N": 14},
                "G": { "-": 9,"A": 13,"T": 13,"C": 13,"G": 6,"N": 14},
                "N": { "-": 14,"A": 14,"T": 14,"C": 14,"G": 14,"N": 14}
                }
        },
    "2": {
        "A" : { "A": 1, "-": 0, "G": 3, "C": 3, "T": 3, "N": 3},
        "T" : { "T": 1, "A": 3, "-": 0, "G": 3, "C": 3, "N": 3},
        "C" : { "C": 1, "T": 3, "A": 3, "-": 0, "G": 3, "N": 3},
        "G" : { "G": 1, "C": 3, "T": 3, "A": 3, "-": 0, "N": 3},
        "-" : { "-": None, "G": 2, "C": 2, "T": 2, "A": 2, "N": 3},
        "N" : { "-": 3, "G": 3, "C": 3, "T": 3, "A": 3, "N": 3}
        }
    }

def observation_probs(index):
    num_obs = 15
    prob = 0.95
    index_prob = prob/len(index)
    other_prob = (1-prob)/(num_obs-len(index)+1)
    probs = [other_prob]*num_obs
    for i in index:
        probs[i] = index_prob
    obs = list(range(0,15))
    return dict(zip(obs,probs))

def get_obs_probs(index,complex_):
    obs_probs = observation_probs([index])
    if complex_:
        obs_probs[index] = 0.94
        obs_probs[0] = 0.01
    return obs_probs

def get_states():
    states = []
    states_with_index = { "INS_1|1": 10,"INS_1|0": 5,"INS_0|1": 9,
                          "DEL_1|1": 12,"DEL_1|0": 6,"DEL_0|1": 7 }
    for state_and_index in states_with_index:
        index = states_with_index[state_and_index]
        for complex_ in [True,False]:
            obs_probs = get_obs_probs(index,complex_)
            name = state_and_index
            if complex_:
                name = "COMPLEX.%s" % name
            state = State(DiscreteDistribution(obs_probs),name=name)
            states.append(state)
    #error_obs_prob = observation_probs([1,2,3,4])  #[4,5,6,7,8,9,11,13])
    #error_state = State(DiscreteDistribution(error_obs_prob),name="ERROR")
    #states.append(error_state)
    complex_obs_prob = observation_probs([5,6,7,9,10,12]) #list(range(0,15)))
    complex_state = State(DiscreteDistribution(complex_obs_prob),name="COMPLEX_.|.")
    states.append(complex_state)
    return states
            
def three_hmm():
    model = HiddenMarkovModel("Sequence Aligner")
    states = get_states()
    normal = State(DiscreteDistribution(observation_probs([0])), name="NORMAL")
    states.append(normal)
    num_states = len(states)
    model.add_states(states)
    for state in states:
        if state.name != "NORMAL":
            model.add_transition(model.start, state, 0)
            model.add_transition(state, model.end, 0)
        else:
            model.add_transition(model.start, state, 1)
            model.add_transition(state, model.end, 1)
    trans_probs = np.zeros((num_states,num_states))
    trans_probs.fill(1e-15)
    np.fill_diagonal(trans_probs,.99999) # Stay in the same state
    trans_probs[:,len(trans_probs)-1] = 1e-15 # Transition to normal
    trans_probs[len(trans_probs)-1] = [4.5e-15]*num_states
    trans_probs[len(trans_probs)-1][-1] = 1-2e-15 # normal event
    #print >> sys.stderr, trans_probs
    for state, t_prob in zip(states,trans_probs):
        for s, prob in zip(states,t_prob):
            model.add_transition(state,s,prob)
            #print >> sys.stderr, state.name,s.name, prob
    model.add_transition(normal,normal,(1-2e-15))
    model.bake()
    return model

model = {"3": three_hmm()}
