cd ..
#for c in vertical horizontal groundmotion
#do
#python plot_stations.py d 2014p240655 "${c}"
#done
for c in vertical horizontal groundmotion                                        
do                                                                               
python plot_stations.py d 2015p822263 "${c}"                                     
done   
#for c in vertical horizontal groundmotion                                        
#do                                                                               
#python plot_stations.py d 2016p892721 "${c}"                                     
#done   
for c in vertical horizontal groundmotion                                        
do                                                                               
python plot_stations.py d 2017p059122 "${c}"                                     
done   
