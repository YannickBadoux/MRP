start_time=$(date +%s)

python ~/Documents/MRP/MRP_code/acs_far_away_distance.py --far_away_dist 30 --output ~/Documents/MRP/acs_far_away_dist_test
python ~/Documents/MRP/MRP_code/acs_far_away_distance.py --far_away_dist 50 --output ~/Documents/MRP/acs_far_away_dist_test
python ~/Documents/MRP/MRP_code/acs_far_away_distance.py --far_away_dist 70 --output ~/Documents/MRP/acs_far_away_dist_test
python ~/Documents/MRP/MRP_code/acs_far_away_distance.py --far_away_dist 100 --output ~/Documents/MRP/acs_far_away_dist_test
python ~/Documents/MRP/MRP_code/acs_far_away_distance.py --far_away_dist 150 --output ~/Documents/MRP/acs_far_away_dist_test

end_time=$(date +%s)
echo "Time elapsed: $((end_time-start_time)) seconds"