import pandas as pd
from tqdm import tqdm_notebook
import cobra
from cobra import Reaction, Metabolite
import operator # for sorting dict
pd.set_option('display.max_rows', 10000) # Show everything
pd.set_option('display.max_colwidth', -1)
pd.set_option('expand_frame_repr', False)

def findBiomarkers(model,fvaRxns=[],mods=[],mode='',metabolomics=False,method=2,synchronous=False,
                eps=1.0,cutoff=0.1,fracOpt=0,forceFlux=True,looseGeneAssociation=False):
    ''' returns a pandas dataframe listing biomarkers in cobrapy model M. under
        alteration of reactions in mods.

    geneIDs:        in IEMgene mode the list of genes to knockout.
    fvaRxns:        FVA is performed on this subset of exchange reactions in M.
    mods:           a list of reactions/genes that will be altered. E.g. the gene/reactions affected by the IEM or the drug import reaction. 
    mode:           we have implemented three methods: IEMrxn (Schlomi et al. method),
                    IEMgene (same with gene names) and drug
    metabolomics    Boolean. If true the metabolites in fvaRxns (change this later to reporters or something) get a temporary sink reaction for which the change in the FVA interval will be reported. 
    method:         indicates whether to take the union of the forward and
                    backward healthy intervals (1) or only take it when the
                    forward interval is [0,0] (2).
    synchronous:    if controlRxns contains multiple entries do everything
                    simultaneously or separately
    eps:            the flux to force through the controlRxns in the healthy case.
                    In the disease state the same reactions are blocked, i.e.
                    carry zero flux.
    cutoff:         the minimal change percentage (value of the Score variable)
                    in the interval size to be considered significant
    fracOpt:        input for FVA, i.e. the minimal percentage of the optimum of
                    the objective function to attain
    forceFlux:      whether or not to force flux through the reaction in the
                    healthy case
    looseGeneAssociation: ?

	The default settings with IEMgene/IEMrxn reproduce the method used by Schlomi et al (2009).  '''

    # interpret input arguments
    # find the mode if it was not given
    if mode == '':
        if all([gene in model.genes for gene in mods]):
            mode = 'IEMgene'
        elif all([rxn in model.reactions for rxn in mods]): # this works even for strings
            mode = 'IEMrxn'
        else:
            if len(mods) == 1:
                try:
                    if model.reactions.get_by_id(mods[0]).products == []:
                        mode = 'drug'
                except Exception as e:
                    raise
            else:
                print('Cannot identify mods input')
                return pd.DataFrame(columns=['ID','Name','Reaction','Prediction',
                                            'WT','Mutant','Score'])

        print('Interpreting mode as {}').format(mode)


    # rename and reinterpret the input based on the mode
    if mode == 'IEMgene':
        geneIDs = mods
        print "Finding affected reactions using gene knockout."
        controlRxns = cobra.manipulation.find_gene_knockout_reactions(model,mods)
        if looseGeneAssociation:
            print "Finding affected reactions using gene-reaction association."
            controlRxns = [rxn for gene in geneIDs for rxn in model.genes.get_by_id(gene).reactions]

        if len(controlRxns) == 0:
            print('This gene list: {}, does not affect any reactions.'.format(geneIDs))
            return pd.DataFrame(columns=['ID','Name','Reaction','Prediction',
                'WT','Mutant','Score'])

    elif mode == 'IEMrxn' or mode == 'drug':
        controlRxns = mods
        if type(controlRxns[0]) == str:
            controlRxns = [ model.reactions.get_by_id(ID) for ID in controlRxns] # take rxn IDs in, turn into rxn objects
    else:
        print 'Cannot parse mode %s. Choose from: IEMgene, IEMrxn and drug.' % (mode)
        return
            
    
    # Are we performing metabolomics or biomarker prediction? 
    # Check if user gave us metabolites and if so add sink reactions to the model
    # set fvaRxns to the list of added sink reactions
    if metabolomics:
        if not(all([met in model.metabolites for met in fvaRxns]) or all([met in [m.id for m in model.metabolites] for met in fvaRxns])):
            print 'For metabolomics mode the fvaRxns input needs to consist of metabolites or their IDs.'
            return
        if fvaRxns == []:
            print 'For metabolomics prediction fvaRxns needs to contain a list of metabolites to simulate with.'
                # print('Empty fvaRxns list given: generating a list of sink/exchange reactions for the most highly connected metabolites in the model.')
                # model,sinks,sink_mets = generateSinkReactions(model,setting='ext')
                # fvaRxns.extend(sinks)
                # fvaRxns = list(set(fvaRxns)) # no duplicates
                # print 'Performing FVA on:'
                # print
        else:
            model,sinks,sink_mets = generateSinkReactions(model,fvaRxns,setting='ext')
            fvaRxns = list(set(sinks)) # no duplicates


    rxnIDs = [rxn.id for rxn in controlRxns ]
    rxnNames = [ rxn.name[0:50] for rxn in controlRxns ]
    rxnReactions = [ rxn.reaction for rxn in controlRxns ]
    df = pd.DataFrame([], index=rxnIDs, columns=['Description','Reaction'])
    df.Description = rxnNames; df.Reaction = rxnReactions

    # print informative message
    print('Modifications will be performed on the following reactions:')
    print df
    print

    # print('FVA will be performed on the following reactions:')
    # print fvaRxns
    # print


    # PERFORM THE ALGORITHM
    # THE LOGIC:
    # 1. calculate the fluxes before and after mutation either reaction by reaction (Schlomi et al.) or all at once
    # 2. combine forward and backward intervals for the healthy case
    # 3. detect biomarkers through changes in WT to mutant intervals
    # 4. generate output dataframe

    # when mutating one reaction at a time the biomarkers need to be consolidated at each step
    if not synchronous and len(controlRxns)>1:
        biomarkerCount = {} # keep track of various (contradicting) predictions
        biomarkerTable = pd.DataFrame(columns=['ID','Name','Reaction','Prediction',
        'WT','Mutant','Score'])
        for rxn in tqdm_notebook(controlRxns):
            [WTf,WTb,mutant] = calcFluxes(model,[rxn],fvaRxns,eps,fracOpt,mode,forceFlux,metabolomics)
            [WTint,mutantint] = uniteForwBack(WTf,WTb,mutant,fvaRxns,method,mode)
            if WTint == {} and mutantint == {}: # nothing to see here
                print('Empty wild-type and mutant FVA intervals. Skipping to next reaction.')
                continue
            [biomarkerRxns,biomarkers,score,extLvl] = predictBiomarkers(model,WTint,mutantint,fvaRxns,cutoff)
            subTable = genTable(biomarkers,biomarkerRxns,score,extLvl,WTint,mutantint,cutoff,synchronous,mode)
            [biomarkerTable,biomarkerCount] = updateTable(subTable,biomarkerTable,biomarkerCount)

    else: # change all controlRxns at once
        [WTf,WTb,mutant] = calcFluxes(model,controlRxns,fvaRxns,eps,fracOpt,mode,forceFlux,metabolomics)
        [WTint,mutantint] = uniteForwBack(WTf,WTb,mutant,fvaRxns,method,mode)
        [biomarkerRxns,biomarkers,score,extLvl] = predictBiomarkers(model,WTint,mutantint,fvaRxns,cutoff)
        biomarkerTable = genTable(biomarkers,biomarkerRxns,score,extLvl,WTint,mutantint,cutoff,synchronous,mode)

    # sort and cut off the biomarker scores
    if cutoff != 0:
        print '{} low confidence biomarkers with scores below the cutoff were found'.format(len(biomarkerTable[biomarkerTable.Score < cutoff]),len(biomarkerTable[biomarkerTable.Score == 0]))
    biomarkerTable = biomarkerTable.sort_values(by='Score',ascending=False)
    significantBiomarkerTable = biomarkerTable[biomarkerTable.Score >= cutoff]

    # formatting
    significantBiomarkerTable = significantBiomarkerTable.reset_index()
    significantBiomarkerTable = significantBiomarkerTable[['ID','Name','Prediction','WT','Mutant','Score']]

    return significantBiomarkerTable

def generateSinkReactions(M,sink_mets,setting):
    def sink_exists(met):
        ''' Find if an exchange for met exists'''
        for r in M.reactions:
            if len(r.metabolites) == 1 and met in r.metabolites:
                return True,r.id
        # if we get here such a reaction was not found
        return False,''  

    # make sure sink_mets are metabolite objects
    if all([type(m) == str or type(m) == unicode for m in sink_mets]):
        sink_mets = [M.metabolites.get_by_id(m) for m in sink_mets]

    # add sink reactions
    sinks = []
    for met in sink_mets:
        b,existing_sink = sink_exists(met)
        if not b:
            rxn = Reaction('tmp_sink_' + met.id)
            M.add_reaction(rxn)
            rxn.add_metabolites({met:-1})
            rxn.upper_bound = 0; rxn.lower_bound = 0 # have to init to zero. Only activate when simulating
            sinks.append(rxn.id)
        else:
            sinks.append(existing_sink)
            print 'Taking existing sink reaction', existing_sink

    return M, sinks, sink_mets

def calcFluxes(model,controlRxns,fvaRxns,eps,fracOpt,mode,forceFlux,metabolomics):
    ''' '''

    def calcFluxes_drug(M,rxnlist):
        # The case with an influx of the control reaction
        # forward
        # WT: without the drain
        for r in controlRxns:
            M.reactions.get_by_id(r.id).lower_bound = 0 # no influx
            M.reactions.get_by_id(r.id).upper_bound = 1000 # possible outflux
        
        WTf = cobra.flux_analysis.variability.flux_variability_analysis(M,reaction_list=rxnlist,fraction_of_optimum=fracOpt)

        # no backward needed, only drain in outward direction

        # mutant: 
        for r in controlRxns:
            M.reactions.get_by_id(r.id).lower_bound = -1000
            M.reactions.get_by_id(r.id).upper_bound = -eps/len(controlRxns)
        
        mutant = cobra.flux_analysis.variability.flux_variability_analysis(M,reaction_list=rxnlist,fraction_of_optimum=fracOpt)

        return WTf,mutant

    def calcFluxes_IEM(M,rxnlist):
        # just run FVA on all rxns at once
        # Healthy case
        # forward
        # M = M_orig.copy() # keep model as is to be able to reset bounds

        if forceFlux:
            for rxn in controlRxns:
                if rxn.upper_bound > 0:
                    M.reactions.get_by_id(rxn.id).lower_bound = eps/len(controlRxns)

        try:
            WTf = cobra.flux_analysis.variability.flux_variability_analysis(M,reaction_list=rxnlist,fraction_of_optimum=fracOpt) # DATAFRAME
        except:
            print('forward FVA could not be solved. Continuing without the forward interval.')
            WTf = pd.DataFrame()

        # backward
        M = model.copy() # reset bounds

        if forceFlux:
            for rxn in controlRxns:
                if rxn.lower_bound < 0:
                    M.reactions.get_by_id(rxn.id).upper_bound = -eps/len(controlRxns)

        try:
            WTb = cobra.flux_analysis.flux_variability_analysis(M,reaction_list=rxnlist,fraction_of_optimum=fracOpt)
        except:
            print('backward FVA could not be solved. Continuing without the backward interval.')
            WTb = pd.DataFrame()

        # Disease case
        M = model.copy() # reset bounds
        for rxn in controlRxns:
            M.reactions.get_by_id(rxn.id).lower_bound = M.reactions.get_by_id(rxn.id).upper_bound = 0
        mutant = cobra.flux_analysis.flux_variability_analysis(M,reaction_list=rxnlist,fraction_of_optimum=fracOpt)

        return WTf,WTb,mutant


    model_orig = model.copy()

    WTf = pd.DataFrame(columns=['minimum','maximum']) # init
    WTb = pd.DataFrame(columns=['minimum','maximum']) # init
    mutant  = pd.DataFrame(columns=['minimum','maximum'])

    if metabolomics: # run fva on single reactions so that we can block the previous ones we added when metabolomics is on
        for rxn in tqdm_notebook(fvaRxns,leave=False):
            with model as model: # every time undo the changes to the bounds and skip using copy for speed

                # open the sink for this metabolite
                model.reactions.get_by_id(rxn).lower_bound = -1000 # no influx
                model.reactions.get_by_id(rxn).upper_bound = 1000 # possible outflux

                if mode == 'drug':
                    subWTf,submutant = calcFluxes_drug(model,[rxn])
                else:
                    subWTf,subWTb,submutant = calcFluxes_IEM(model,[rxn])

                # update the total set
                WTf.loc[rxn] = subWTf.loc[rxn]
                mutant.loc[rxn] = submutant.loc[rxn]
                if mode != 'drug':
                    WTb.loc[rxn] = subWTb.loc[rxn]
            
    else: # we can run all fvaRxns at once no need to change bounds in the middle
        with model as model:
            if mode == 'drug':
                WTf,mutant = calcFluxes_drug(model,fvaRxns)
            else:
                WTf,WTb,mutant = calcFluxes_IEM(model,fvaRxns)
        
    return WTf, WTb, mutant

def uniteForwBack(WTf,WTb,mutant,fvaRxns,method,mode):
    ''' '''

    # take union of forward and backward healthy intervals
    WT = pd.DataFrame(columns=WTf.columns) # init

    if len(WTf) != 0 and len(WTb) == 0: # no backward calculation
        WT = WTf
    elif len(WTb) != 0 and len(WTf) == 0: # no forward calculation
        WT = WTb
    elif len(WTb) == 0 and len(WTf) == 0: # no results
        print('No healthy interval could be calculated')
    elif method == 1 and len(WTf) != 0 and len(WTb) != 0:
        for rxn in fvaRxns: # enlarge interval with WTb if needed
            if WTb.loc[rxn]["minimum"] < WTf.loc[rxn]["minimum"]:
                WT.loc[rxn]["minimum"] = WTb.loc[rxn]["minimum"]    
            if WTb.loc[rxn]["maximum"] > WTf.loc[rxn]["maximum"]:
                WT.loc[rxn]["maximum"] = WTb.loc[rxn]["maximum"]
    elif method == 2 and len(WTf) != 0 and len(WTb) != 0:
        for rxn in fvaRxns: # enlarge interval with WTb if needed
            if WTf.loc[rxn]["minimum"] == WTf.loc[rxn]["maximum"] == 0:
                WT.loc[rxn] = WTb.loc[rxn]
            else:
                WT.loc[rxn] = WTf.loc[rxn]
    else:
        print('Something weird is going on!')
        return

    # analyse the results and predict biomarkers
    WTint = {}
    mutantint = {}

    if len(WT) != 0 and len(mutant) != 0:
        for r in fvaRxns:
            WTint[r] = [ round(WT.loc[r]["minimum"],3), round(WT.loc[r]["maximum"],3) ]
            mutantint[r] = [ round(mutant.loc[r]["minimum"],3), round(mutant.loc[r]["maximum"],3) ]

    return WTint, mutantint

def predictBiomarkers(M,WTint,mutantint,fvaRxns,cutoff):
    ''' '''

    # be very careful when interpreting these ranges: it is WT - mutant.
    # A negative lb difference means the mutant has reduced uptake capabilities which results in higher serum levels.
    # Positive lb difference occurs when mutants have increased uptake capabilities meaning lower serum levels.
    # positive ub difference means mutant can produce less, so lower serum levels.
    # negative ub difference means mutants can produce more so higher serum levels.
    score = {}
    extLvl = {}
    for rxn in fvaRxns:
        # calculate score
        lb = [WTint[rxn][0],mutantint[rxn][0]]; ub = [WTint[rxn][1],mutantint[rxn][1]]

        if lb == [0,0]:
            change_lower_bound = 0
        else:
            change_lower_bound = abs(max(lb)-min(lb))/max([abs(el) for el in lb])
        if ub == [0,0]:
            change_upper_bound = 0
        else:
            change_upper_bound = abs(max(ub)-min(ub))/max([abs(el) for el in ub])

        score[fvaRxns.index(rxn)] = max( change_lower_bound,change_upper_bound )

        # determine direction of change
        if WTint[rxn][0] == mutantint[rxn][0] and WTint[rxn][1] == mutantint[rxn][1]:
            extLvl[rxn] = "Unchanged"
        elif WTint[rxn][1] < mutantint[rxn][0]:
            extLvl[rxn] = "H.C. Elevated"
        elif WTint[rxn][0] > mutantint[rxn][1]:
            extLvl[rxn] = "H.C. Reduced"
        elif WTint[rxn][0] <= mutantint[rxn][0] and WTint[rxn][1] <= mutantint[rxn][1] and ( max( abs(WTint[rxn][0] - mutantint[rxn][0]), abs(WTint[rxn][1] - mutantint[rxn][1]) ) > 0 ) :
            extLvl[rxn] = "Elevated"
        elif WTint[rxn][0] >= mutantint[rxn][0] and WTint[rxn][1] >= mutantint[rxn][1] and ( abs(WTint[rxn][0] - mutantint[rxn][0]) > 0 or abs(WTint[rxn][1] - mutantint[rxn][1]) > 0 ) :
            extLvl[rxn] = "Reduced"
        else:
            extLvl[rxn] = "Undetermined"

    # identify changed intervals indicating possible biomarkers
    if cutoff > 0:
        biomarkerRxns = [ M.reactions.get_by_id(rxn) for rxn in fvaRxns if extLvl[rxn] !=  'Unchanged' ]
    else:
        biomarkerRxns = [ M.reactions.get_by_id(rxn) for rxn in fvaRxns]
    score = [ score[fvaRxns.index(rxn.id)] for rxn in biomarkerRxns ]
    biomarkers = [ rxn.metabolites.keys() for rxn in biomarkerRxns ] # this is a list of lists
    biomarkers = [ item for sublist in biomarkers for item in sublist ]

    return biomarkerRxns, biomarkers, score, extLvl

def updateTable(subTable,biomarkerTable,biomarkerCount):
    '''Fix duplicates. Check if biomarker already exists.
    If so, check if the qualitative prediction is the same.
    If it is, keep it, if it is not delete the biomarker.
    This leads to a majority rule scenario '''
    for bm in subTable['ID'].tolist():
        currentRow = subTable.loc[subTable['ID'] == bm]

        # update the count
        if bm not in biomarkerCount.keys(): biomarkerCount[bm] = 0
        if 'Elevated' in currentRow['Prediction'].tolist() or 'H.C. Elevated' in currentRow['Prediction'].tolist():
            biomarkerCount[bm] += 1
        elif 'Reduced' in currentRow['Prediction'].tolist() or 'H.C. Reduced' in currentRow['Prediction'].tolist():
            biomarkerCount[bm] -= 1

        # based on the count choose to keep or drop the biomarker
        if biomarkerCount[bm] == 0 and bm in biomarkerTable['ID'].tolist(): # equal contradictory predictions
            print('Removed {} because it has an equal number of contradictory predictions.'.format(bm))
            biomarkerTable = biomarkerTable[biomarkerTable['ID'] != bm]
        elif biomarkerCount[bm] != 0 and bm not in biomarkerTable['ID'].tolist():
            biomarkerTable = biomarkerTable.append(currentRow) # add new biomarker

    return biomarkerTable,biomarkerCount

def genTable(biomarkers,biomarkerRxns,score,extLvl,WTint,mutantint,cutoff,synchronous,mode):
    ''' '''

    # if mode == 'metabolomics':
    #     # flip WT to be without drain and mutant to be with drain
    #     # also flip the predictions
    #     tmp = WTint
    #     WTint = mutantint; mutantint = tmp 

    #     flip = {'Elevated':'Reduced','Reduced':'Elevated','H.C. Elevated':'H.C. Reduced','H.C. Reduced':'H.C. Elevated',
    #     'Unchanged':'Unchanged','Undetermined':'Undetermined'}
    #     extLvl = {r:flip[extLvl[r]] for r in extLvl}


    # Generate the output pandas dataframe, either return or print it
    biomarkerTable = pd.DataFrame({'ID': [ bm.id for bm in biomarkers ],
                'Name': [ bm.name for bm in biomarkers ],
                'Reaction': [ rxn.id for rxn in biomarkerRxns ],
                'Prediction': [ extLvl[rxn.id] for rxn in biomarkerRxns ],
                'WT': [ WTint[rxn.id] for rxn in biomarkerRxns ],
                'Mutant': [ mutantint[rxn.id] for rxn in biomarkerRxns ],
                'Score': score })

    return biomarkerTable[['ID','Name','Reaction','Prediction','WT','Mutant','Score']].round(3)
# end
