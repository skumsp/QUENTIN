function DSamp_dir = findDSampDir1(DSamp)
AMSamp_dir = (DSamp <= DSamp');
DSamp_dir = DSamp.*AMSamp_dir;