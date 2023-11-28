mas_up = np.copy(tmplate)
spots  = 10
span   = 30

# original guesses for peak time offsets
offs_peak = [[-7,0,7],[-14,-7,0,7,14],[-7,0,7]]

# we want to start from about Nov23
# loop over masers
M = []
pmax = mas_up.shape[0]

mas_array = np.arange(len(masers));

for shuffles in range(10):
    # start with zero fails
    fail = 0
    #random start maser and order
    np.random.shuffle(y); 
    # zero results
    mas_up = np.copy(tmplate)
    mpeak  = np.zeros(shape=(len(masers),11))
    # loop over masers from shuffled array
    for mas in mas_array:
        ra  = masers[mas][0]
        # have a few attempts per maser
        for mas_trys in range(10):
            # loop over 3 seasons per maser, start at epoch 1
            ep = 0
            for season in range(3):     
                # get ideal peaks from file
                pyr,mnt,day = masers[mas][1:-1][season]
                # loop over epochs in offs_peak
                n = 0
                for k in offs_peak[season]:
                    n += 1
                    # for each epoch, try 1000 spots
                    for ep_try in range(spots):
                        # calculate the peak day in UT since start of 2023
                        pday = mnth_day[mnt-1+12*pyr]+day+k
                        # making sure they're after 24th Nov 2023, else fail
                        if pday<=335:
                            # new spot position
                            k = int(np.random.uniform(-span,span))
                            continue
                        # making sure they're in the 3 year range, else fail
                        if pday>=pmax:     
                            # new spot position
                            k = int(np.random.uniform(-span,span))
                            continue

                        #check if the day/time is reserved by IVS
                        ivs_res = np.any(year_res[pday,(abs(sid_ced[pday%365,:]-ra)<4)].astype(bool)) 
                        # check if the day is already reserved by spirals
                        spi_res = np.any(mas_up[pday,(abs(sid_ced[pday%365,:]-ra)<4)].astype(bool))            
                        if ~ivs_res and ~spi_res:
                            #ep += 1
                            mas_up[pday,(abs(sid_ced[pday%365,:]-ra)<4)] = 3+i
                            mpeak[mas,ep] = pday
                            break 
                        else: 
                            # new spot position
                            #if ivs_res and spi_res: res = 'both:'
                            #if ivs_res and ~spi_res: res = 'ivs:'                        
                            #if ~ivs_res and spi_res: res = 'spi:'                        
                            k = int(np.random.uniform(-span,span))
                            continue
                    #if we get to here we've either got a solution or failed to converge
                    # if we've failed a spot, might as well stop
                    if ep_try==spots-1: 
                        pday = mnth_day[mnt-1+12*pyr]+day+offs_peak[j][n-1]
                        mpeak[mas,ep] = pday
                        mas_up[pday,(abs(sid_ced[pday%365,:]-ra)<4)] = 3+i
                        fail += 1
                    ep += 1

