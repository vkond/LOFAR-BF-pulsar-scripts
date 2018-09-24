h5dump -a OBSERVATION_START_UTC $@ | grep "(0)" | sed 's|   (0): \"||g' | sed 's:.000000000Z\"::g' | sed 's:T: :g'
