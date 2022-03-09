#include <ilcplex/ilocplex.h> 
#include <stdio.h>   
#include <stdlib.h>            
#include <sys/time.h>              
#include <stdlib.h>            
#include <time.h>                         
#include <fstream>  
#include <string>
#include <math.h>
#include <sstream>
#include <algorithm>           
#include <vector>
#include <omp.h> 
#include "Avgminmax02.h"
#include "Seqinsertion.h"
#define random(x) (rand()%x)
namespace patch {template < typename T > std::string to_string( const T& n ){std::ostringstream stm; stm << n; return stm.str() ;}}

using namespace std;
//input data starts here
extern const int CNN = 1;
extern const int P = 3;
extern const int PT = 6;
extern const int T = 18;
extern const int I1 =1;
extern const int I2 = 8;
extern const int I = 9;
extern const int BN = 72;
extern const int LBN = 9;
extern const int SPN = 340;
extern const int SN = 10;
extern const int MBSN = 9;
double PRO[P][SN]={{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1},{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1},{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1}}; 
int PTN[P]={6,6,6}; 
int PTset[P][PT]={{0,1,2,3,4,5},{6,7,8,9,10,11},{12,13,14,15,16,17}}; 
int sup[I][P][SN][PT]={{{{11589,0,0,0,0,0},{11615,0,0,0,0,0},{11589,0,0,0,0,0},{11623,0,0,0,0,0},{11604,0,0,0,0,0},{11495,0,0,0,0,0},{11495,0,0,0,0,0},{11561,0,0,0,0,0},{11710,0,0,0,0,0},{11608,0,0,0,0,0}},{{11524,0,0,0,0,0},{11786,0,0,0,0,0},{11659,0,0,0,0,0},{11622,0,0,0,0,0},{11506,0,0,0,0,0},{11518,0,0,0,0,0},{11617,0,0,0,0,0},{11380,0,0,0,0,0},{11594,0,0,0,0,0},{11473,0,0,0,0,0}},{{11663,0,0,0,0,0},{11354,0,0,0,0,0},{11759,0,0,0,0,0},{11577,0,0,0,0,0},{11407,0,0,0,0,0},{11519,0,0,0,0,0},{11529,0,0,0,0,0},{11512,0,0,0,0,0},{11619,0,0,0,0,0},{11535,0,0,0,0,0}}},{{{-217,-185,-201,-209,-205,-213},{-221,-215,-219,-187,-219,-217},{-215,-215,-191,-187,-217,-207},{-193,-221,-221,-201,-189,-203},{-185,-215,-189,-207,-185,-213},{-219,-201,-223,-215,-189,-193},{-201,-199,-207,-199,-191,-209},{-189,-189,-199,-201,-207,-207},{-207,-191,-211,-201,-215,-187},{-219,-207,-207,-209,-201,-207}},{{-187,-207,-197,-197,-219,-185},{-201,-195,-189,-215,-219,-219},{-205,-187,-203,-223,-189,-205},{-221,-217,-215,-221,-205,-201},{-189,-221,-187,-193,-219,-199},{-193,-211,-221,-191,-211,-191},{-217,-195,-217,-215,-199,-211},{-201,-215,-199,-189,-183,-201},{-199,-211,-199,-213,-209,-191},{-215,-187,-223,-215,-183,-221}},{{-215,-205,-191,-217,-197,-189},{-183,-207,-199,-191,-217,-205},{-201,-183,-185,-197,-219,-221},{-213,-189,-199,-199,-199,-219},{-207,-183,-187,-205,-215,-221},{-191,-199,-183,-187,-209,-199},{-213,-199,-197,-185,-193,-215},{-197,-199,-199,-187,-201,-219},{-219,-183,-199,-223,-203,-193},{-185,-221,-187,-207,-215,-221}}},{{{-171,-175,-158,-163,-147,-155},{-172,-145,-156,-155,-171,-150},{-159,-151,-153,-169,-153,-153},{-153,-155,-169,-177,-158,-156},{-174,-161,-155,-172,-166,-177},{-172,-174,-169,-175,-155,-169},{-172,-150,-147,-148,-171,-158},{-158,-147,-164,-159,-156,-151},{-175,-172,-150,-150,-172,-169},{-145,-150,-171,-145,-159,-163}},{{-161,-153,-163,-155,-155,-151},{-172,-172,-161,-164,-150,-174},{-174,-175,-172,-150,-156,-174},{-171,-161,-172,-169,-150,-161},{-147,-148,-161,-159,-171,-167},{-177,-145,-151,-166,-155,-169},{-175,-161,-159,-177,-148,-167},{-153,-169,-147,-159,-151,-163},{-156,-172,-151,-172,-153,-174},{-151,-174,-166,-164,-175,-151}},{{-148,-151,-166,-151,-172,-167},{-167,-148,-174,-169,-158,-174},{-169,-177,-163,-151,-177,-166},{-174,-151,-171,-177,-148,-150},{-158,-158,-161,-167,-177,-166},{-174,-177,-148,-161,-151,-156},{-147,-175,-148,-172,-177,-148},{-155,-155,-147,-166,-158,-175},{-175,-167,-172,-155,-147,-166},{-171,-174,-161,-153,-164,-147}}},{{{-229,-219,-249,-224,-222,-233},{-249,-215,-226,-236,-243,-254},{-210,-236,-252,-224,-215,-236},{-249,-231,-217,-243,-236,-254},{-240,-252,-226,-238,-240,-231},{-245,-249,-247,-226,-217,-226},{-249,-231,-233,-213,-208,-219},{-245,-217,-236,-252,-215,-215},{-243,-238,-238,-249,-243,-231},{-217,-231,-245,-245,-215,-254}},{{-217,-243,-252,-217,-240,-233},{-226,-213,-245,-254,-238,-224},{-254,-243,-215,-247,-215,-233},{-238,-243,-222,-236,-243,-217},{-247,-224,-254,-215,-224,-226},{-249,-229,-226,-247,-245,-254},{-252,-247,-247,-240,-226,-240},{-222,-226,-213,-219,-217,-224},{-208,-247,-252,-231,-240,-252},{-208,-210,-210,-254,-245,-224}},{{-254,-226,-245,-226,-247,-229},{-208,-222,-245,-247,-222,-245},{-247,-236,-219,-245,-226,-233},{-222,-217,-233,-247,-236,-222},{-238,-229,-219,-245,-215,-215},{-213,-252,-208,-208,-249,-249},{-249,-254,-208,-243,-226,-215},{-219,-213,-231,-245,-213,-215},{-231,-243,-245,-219,-243,-238},{-217,-231,-219,-249,-226,-245}}},{{{-224,-215,-197,-206,-202,-210},{-215,-230,-221,-215,-202,-239},{-224,-206,-215,-232,-202,-200},{-195,-215,-200,-239,-221,-232},{-226,-208,-237,-195,-221,-195},{-197,-228,-195,-230,-202,-208},{-239,-213,-206,-228,-234,-219},{-195,-230,-221,-228,-237,-219},{-234,-226,-197,-200,-224,-206},{-195,-213,-217,-210,-202,-224}},{{-204,-219,-237,-234,-208,-204},{-228,-226,-215,-221,-234,-226},{-232,-237,-213,-221,-195,-232},{-206,-217,-226,-226,-195,-228},{-195,-237,-217,-239,-210,-217},{-237,-204,-204,-206,-239,-217},{-210,-202,-221,-224,-224,-215},{-234,-228,-224,-217,-224,-234},{-204,-217,-230,-210,-234,-239},{-200,-202,-202,-226,-224,-230}},{{-200,-197,-234,-230,-224,-215},{-215,-230,-230,-234,-200,-224},{-237,-232,-224,-204,-237,-234},{-228,-200,-228,-202,-221,-200},{-232,-204,-219,-228,-237,-226},{-213,-197,-217,-210,-230,-217},{-226,-206,-206,-232,-230,-206},{-197,-208,-206,-210,-204,-202},{-239,-234,-213,-208,-197,-237},{-224,-228,-200,-208,-219,-200}}},{{{-269,-253,-287,-285,-258,-266},{-261,-261,-255,-261,-255,-282},{-285,-266,-266,-263,-271,-285},{-279,-255,-277,-287,-290,-277},{-255,-258,-290,-245,-290,-274},{-239,-269,-239,-245,-269,-279},{-269,-271,-266,-250,-293,-274},{-277,-242,-247,-277,-290,-274},{-290,-266,-282,-290,-250,-282},{-282,-245,-242,-253,-285,-287}},{{-285,-261,-250,-285,-285,-253},{-258,-290,-290,-250,-266,-258},{-266,-277,-245,-287,-277,-250},{-250,-255,-277,-269,-247,-247},{-277,-266,-266,-282,-287,-287},{-263,-274,-271,-263,-274,-266},{-261,-274,-263,-290,-242,-285},{-239,-282,-250,-247,-266,-274},{-293,-274,-255,-261,-258,-274},{-245,-274,-290,-287,-242,-239}},{{-293,-279,-279,-287,-247,-253},{-245,-239,-242,-253,-290,-245},{-285,-266,-258,-293,-287,-285},{-274,-282,-287,-245,-266,-263},{-253,-247,-242,-261,-239,-293},{-261,-269,-242,-277,-279,-293},{-274,-282,-271,-290,-253,-255},{-261,-266,-269,-293,-287,-250},{-287,-263,-277,-282,-239,-282},{-290,-245,-245,-245,-271,-293}}},{{{-252,-294,-286,-283,-305,-297},{-283,-291,-302,-272,-286,-291},{-280,-308,-291,-286,-294,-305},{-274,-263,-255,-308,-300,-252},{-300,-269,-263,-280,-291,-294},{-291,-277,-297,-302,-291,-255},{-263,-291,-291,-305,-266,-263},{-308,-263,-252,-252,-291,-300},{-300,-291,-305,-308,-269,-263},{-308,-274,-266,-272,-260,-297}},{{-294,-258,-280,-266,-274,-266},{-266,-286,-302,-308,-274,-302},{-283,-255,-305,-291,-266,-274},{-277,-283,-294,-308,-266,-297},{-266,-263,-297,-291,-283,-269},{-260,-274,-277,-255,-266,-269},{-286,-255,-263,-274,-291,-263},{-294,-269,-255,-291,-286,-272},{-308,-260,-302,-266,-255,-255},{-286,-288,-252,-300,-266,-280}},{{-308,-291,-305,-252,-277,-258},{-266,-258,-258,-266,-302,-252},{-302,-297,-283,-308,-269,-252},{-297,-288,-291,-272,-294,-274},{-294,-260,-297,-291,-263,-269},{-252,-258,-294,-300,-305,-297},{-297,-252,-291,-252,-255,-263},{-308,-283,-300,-294,-294,-260},{-288,-294,-302,-274,-258,-269},{-280,-280,-291,-286,-291,-258}}},{{{-302,-300,-300,-272,-269,-288},{-266,-255,-280,-266,-305,-302},{-288,-252,-260,-300,-277,-283},{-272,-294,-258,-300,-288,-297},{-286,-305,-300,-286,-286,-269},{-283,-274,-274,-255,-300,-263},{-272,-255,-305,-294,-288,-277},{-300,-305,-258,-302,-308,-263},{-263,-280,-308,-277,-280,-280},{-294,-277,-305,-266,-305,-305}},{{-266,-294,-305,-272,-300,-272},{-252,-308,-286,-305,-286,-274},{-305,-255,-288,-302,-258,-255},{-308,-286,-280,-291,-272,-297},{-274,-280,-266,-274,-263,-288},{-294,-286,-274,-258,-280,-288},{-269,-258,-255,-272,-294,-263},{-288,-272,-297,-300,-260,-280},{-274,-266,-300,-280,-252,-274},{-300,-263,-294,-283,-280,-266}},{{-286,-286,-297,-302,-274,-297},{-260,-286,-283,-291,-280,-255},{-263,-266,-291,-266,-300,-294},{-277,-294,-308,-305,-255,-305},{-283,-260,-266,-280,-297,-255},{-308,-266,-280,-269,-272,-302},{-294,-260,-305,-269,-297,-308},{-260,-283,-300,-255,-294,-263},{-263,-308,-291,-277,-308,-252},{-272,-305,-255,-266,-305,-263}}},{{{-305,-263,-308,-255,-305,-258},{-291,-294,-283,-255,-269,-277},{-274,-297,-274,-297,-263,-302},{-274,-294,-283,-255,-300,-263},{-305,-272,-252,-258,-274,-294},{-291,-258,-269,-291,-283,-277},{-274,-263,-283,-300,-283,-255},{-277,-291,-308,-302,-277,-260},{-288,-297,-263,-277,-272,-300},{-297,-308,-288,-258,-300,-283}},{{-283,-308,-255,-294,-297,-283},{-294,-277,-294,-308,-277,-294},{-272,-305,-300,-291,-302,-300},{-252,-305,-277,-277,-277,-269},{-280,-269,-272,-263,-294,-283},{-302,-266,-269,-274,-286,-291},{-297,-302,-291,-286,-308,-286},{-260,-272,-266,-297,-291,-260},{-260,-305,-263,-302,-305,-288},{-269,-280,-300,-263,-308,-258}},{{-297,-255,-277,-302,-288,-277},{-280,-269,-297,-274,-277,-272},{-305,-283,-294,-280,-286,-263},{-291,-283,-252,-269,-308,-252},{-272,-266,-300,-272,-283,-252},{-294,-294,-280,-288,-258,-283},{-260,-280,-286,-308,-297,-260},{-302,-297,-302,-274,-308,-288},{-305,-255,-286,-258,-252,-300},{-308,-283,-288,-269,-305,-269}}}}; 
int iniv[I]={2070,261,23,150,170,304,400,300,160}; 
int ubiv[I]={2740,522,414,300,341,608,800,600,320}; 
int len[I1][I2]={{3,6,2,2,3,3,2,1}}; 
double c1[I]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; 
double c2[I]={0.05,7.19,13.32,4.32,2.64,5.16,6.85,4.63,1.82}; 
double c3[BN]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; 
double c4[I1][I2]={{6.54,12.11,3.93,2.4,4.69,6.23,4.21,1.65}}; 
int RLBN[I1][I2]={{9,9,9,9,9,9,9,9}}; 
int LBset[I1][I2][LBN]={{{0,1,2,3,4,5,6,7,8},{9,10,11,12,13,14,15,16,17},{18,19,20,21,22,23,24,25,26},{27,28,29,30,31,32,33,34,35},{36,37,38,39,40,41,42,43,44},{45,46,47,48,49,50,51,52,53},{54,55,56,57,58,59,60,61,62},{63,64,65,66,67,68,69,70,71}}}; 
double frt[BN]={5.23,5.23,5.23,4.58,4.58,4.58,3.92,3.92,3.92,9.69,9.69,9.69,8.48,8.48,8.48,7.27,7.27,7.27,3.14,3.14,3.14,2.75,2.75,2.75,2.36,2.36,2.36,1.92,1.92,1.92,1.68,1.68,1.68,1.44,1.44,1.44,3.75,3.75,3.75,3.28,3.28,3.28,2.81,2.81,2.81,4.98,4.98,4.98,4.36,4.36,4.36,3.74,3.74,3.74,3.37,3.37,3.37,2.95,2.95,2.95,2.53,2.53,2.53,1.32,1.32,1.32,1.16,1.16,1.16,0.99,0.99,0.99}; 
int lbcap[BN]={75,75,75,1,151,151,151,226,226,75,75,75,301,151,151,151,226,226,75,75,75,301,151,151,151,226,226,75,75,75,301,151,151,151,226,226,75,75,75,301,151,151,151,226,226,75,75,75,301,151,151,151,226,226,75,75,75,301,151,151,151,226,226,75,75,75,301,151,151,151,226,226}; 
int ubcap[BN]={150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300}; 
int SHN[BN]={7,4,3,7,4,3,7,4,3,6,3,2,6,3,2,6,3,2,8,4,3,8,4,3,7,4,3,8,4,3,7,4,3,8,4,3,7,4,3,8,4,3,7,4,3,7,4,3,7,4,3,7,4,3,8,4,3,8,4,3,8,4,3,9,4,3,8,4,3,8,5,3}; 
int SHset[BN][MBSN]={{0,1,2,3,4,5,6},{7,8,9,10},{11,12,13},{14,15,16,17,18,19,20},{21,22,23,24},{25,26,27},{28,29,30,31,32,33,34},{35,36,37,38},{39,40,41},{42,43,44,45,46,47},{48,49,50},{51,52},{53,54,55,56,57,58},{59,60,61},{62,63},{64,65,66,67,68,69},{70,71,72},{73,74},{75,76,77,78,79,80,81,82},{83,84,85,86},{87,88,89},{90,91,92,93,94,95,96,97},{98,99,100,101},{102,103,104},{105,106,107,108,109,110,111},{112,113,114,115},{116,117,118},{119,120,121,122,123,124,125,126},{127,128,129,130},{131,132,133},{134,135,136,137,138,139,140},{141,142,143,144},{145,146,147},{148,149,150,151,152,153,154,155},{156,157,158,159},{160,161,162},{163,164,165,166,167,168,169},{170,171,172,173},{174,175,176},{177,178,179,180,181,182,183,184},{185,186,187,188},{189,190,191},{192,193,194,195,196,197,198},{199,200,201,202},{203,204,205},{206,207,208,209,210,211,212},{213,214,215,216},{217,218,219},{220,221,222,223,224,225,226},{227,228,229,230},{231,232,233},{234,235,236,237,238,239,240},{241,242,243,244},{245,246,247},{248,249,250,251,252,253,254,255},{256,257,258,259},{260,261,262},{263,264,265,266,267,268,269,270},{271,272,273,274},{275,276,277},{278,279,280,281,282,283,284,285},{286,287,288,289},{290,291,292},{293,294,295,296,297,298,299,300,301},{302,303,304,305},{306,307,308},{309,310,311,312,313,314,315,316},{317,318,319,320},{321,322,323},{324,325,326,327,328,329,330,331},{332,333,334,335,336},{337,338,339}}; 
int SHsts[SPN]={2,4,6,8,10,12,14,1,5,9,13,0,6,12,1,3,5,7,9,11,13,1,5,9,13,0,6,12,1,3,5,7,9,11,13,0,4,8,12,2,8,14,0,2,4,6,8,10,1,5,9,0,6,0,2,4,6,8,10,1,5,9,0,6,0,2,4,6,8,10,2,6,10,2,8,0,2,4,6,8,10,12,14,2,6,10,14,2,8,14,1,3,5,7,9,11,13,15,1,5,9,13,1,7,13,2,4,6,8,10,12,14,2,6,10,14,1,7,13,1,3,5,7,9,11,13,15,1,5,9,13,2,8,14,2,4,6,8,10,12,14,2,6,10,14,1,7,13,1,3,5,7,9,11,13,15,0,4,8,12,1,7,13,2,4,6,8,10,12,14,2,6,10,14,0,6,12,0,2,4,6,8,10,12,14,0,4,8,12,1,7,13,1,3,5,7,9,11,13,0,4,8,12,1,7,13,2,4,6,8,10,12,14,2,6,10,14,1,7,13,2,4,6,8,10,12,14,2,6,10,14,1,7,13,2,4,6,8,10,12,14,2,6,10,14,0,6,12,1,3,5,7,9,11,13,15,1,5,9,13,0,6,12,1,3,5,7,9,11,13,15,0,4,8,12,1,7,13,0,2,4,6,8,10,12,14,1,5,9,13,1,7,13,0,2,4,6,8,10,12,14,16,1,5,9,13,1,7,13,1,3,5,7,9,11,13,15,1,5,9,13,2,8,14,1,3,5,7,9,11,13,15,0,4,8,12,16,2,8,14}; 
int SHets[SPN]={5,7,9,11,13,15,17,4,8,12,16,3,9,15,4,6,8,10,12,14,16,4,8,12,16,3,9,15,4,6,8,10,12,14,16,3,7,11,15,5,11,17,6,8,10,12,14,16,7,11,15,6,12,6,8,10,12,14,16,7,11,15,6,12,6,8,10,12,14,16,8,12,16,8,14,2,4,6,8,10,12,14,16,4,8,12,16,4,10,16,3,5,7,9,11,13,15,17,3,7,11,15,3,9,15,4,6,8,10,12,14,16,4,8,12,16,3,9,15,3,5,7,9,11,13,15,17,3,7,11,15,4,10,16,4,6,8,10,12,14,16,4,8,12,16,3,9,15,3,5,7,9,11,13,15,17,2,6,10,14,3,9,15,5,7,9,11,13,15,17,5,9,13,17,3,9,15,3,5,7,9,11,13,15,17,3,7,11,15,4,10,16,4,6,8,10,12,14,16,3,7,11,15,4,10,16,5,7,9,11,13,15,17,5,9,13,17,4,10,16,5,7,9,11,13,15,17,5,9,13,17,4,10,16,5,7,9,11,13,15,17,5,9,13,17,3,9,15,3,5,7,9,11,13,15,17,3,7,11,15,2,8,14,3,5,7,9,11,13,15,17,2,6,10,14,3,9,15,2,4,6,8,10,12,14,16,3,7,11,15,3,9,15,1,3,5,7,9,11,13,15,17,2,6,10,14,2,8,14,2,4,6,8,10,12,14,16,2,6,10,14,3,9,15,2,4,6,8,10,12,14,16,1,5,9,13,17,3,9,15}; 

//input data ends here

double t_p4;
int ESTV=10000; 

int TPH[T]; 

extern const int ND1=I1*T; 
extern const int ND2=I2*T; 
extern const int ND=I*T; 
int Nodeset[2][ND]; 

int NDindex[I][T]; 
int PNN[P];  
vector<vector<int>> PNset(P); 
int PNDindex[I][PT]; 
int WND[P]; 
int WNDindex[I2][PT]; 
double AVGD[ND];  

double COST1[ND]; 
double COST2[ND]; 
vector<vector<double>> D;  

double TQ[ND]; 
double IQ[ND]; 

int ARC=0; 
vector<vector<int>> ARCND(2); 
vector<double> COST3; 

int PAN[P];  
vector<vector<int>> PAset(P); 
int PBAN[P][BN];  
vector<vector<vector<int>>> PBAset(P); 
int PINAN[P][ND];
vector<vector<vector<int>>> PINAset(P); 
int POUTAN[P][ND];
vector<vector<vector<int>>> POUTAset(P); 


int PSAN[P]; 
vector<vector<int>> PSAset(P); 
vector<vector<int>> PSAindex(P); 


int INAN[ND]; 
vector<vector<int>> INAset(ND); 
int OUTAN[ND]; 
vector<vector<int>> OUTAset(ND); 


double mfb[BN]; 
vector<int> ARCPN;  
vector<int> ARCPS;  


int titn=0; 
double bidcost;
double thebidlb; 
double transcost;
double stgap=100; 

double Bcap[BN];    
vector<double> Arcflow; 


int SPANN[P];                    
vector<vector<int>> SPANset(P-1);  
vector<vector<double>> MINIdmd(P-1);   
int SPmark[P][ND];     


extern const int TTRD=16; 
int NSOL[TTRD];
vector<vector<double>> VAlsol(TTRD);
vector<vector<double>> PD(TTRD);  
vector<vector<vector<int>>> RSCE(TTRD);   
vector<vector<vector<double>>> PInven(TTRD);   
vector<vector<vector<double>>> PUVval(TTRD);   
vector<vector<vector<double>>> PArcflow(TTRD); 


int SOARCN[P]; 
int SONDN[P];  
vector<vector<int>> SOARCset(P); 
vector<vector<int>> SONDset(P);   
vector<vector<int>> SOARCind (P); 
vector<vector<int>> SONDind (P); 

int MUT=0; 
double LL=9999999;
double UB;
double LB;
double totalgap;
int ALLSN;
double upperbound;

double TTspotcost=0;  
double TTivcost=0;    
double TTbgcost=0;    
double TTbidcap=0;   
double TTbidvol=0;    
double TTspotvol=0;   
double TTivvol=0;    
double TTbgvol=0;    
double TTbirate=0;   

double STspotcost[TTRD];  
double STivcost[TTRD];    
double STbgcost[TTRD];    
double STbidvol[TTRD];    
double STspotvol[TTRD];   
double STivvol[TTRD];    
double STbgvol[TTRD];    
double STbirate[TTRD];   

string num2str(double i)
{    stringstream ss;
ss<<i;
return ss.str();}

void Netsetup(){ 
	for (int p=0;p!=P;++p){
		for (int n=0;n!=PTN[p];++n){
			int thetime=PTset[p][n];
			TPH[thetime]=p;   
		}
	}
	for (int p=0;p!=P;++p){
		PNN[p]=0;
	}
	for (int pp=0;pp!=P-1;++pp){
		SPANN[pp]=0;    
		for (int nd=0;nd!=ND;++nd){
			SPmark[pp][nd]=-1;}
	}
	int counter1=0; 
	for (int i=0;i!=I1;++i){ 
		for (int t=0;t!=T;++t){
			Nodeset[0][counter1]=i;
			Nodeset[1][counter1]=t;
			COST1[counter1]=c1[i]; 
			COST2[counter1]=c2[i]; 
			TQ[counter1]=ubiv[i];  
			IQ[counter1]=iniv[i];  
			NDindex[i][t]=counter1;
			counter1+=1;
		}
	}
	for (int i=0;i!=I2;++i){ 
		for (int t=0;t!=T;++t){
			Nodeset[0][counter1]=i+I1;
			Nodeset[1][counter1]=t;
			COST1[counter1]=c1[i+I1]; 
			COST2[counter1]=c2[i+I1]; 
			TQ[counter1]=ubiv[i+I1];  
			IQ[counter1]=iniv[i+I1];  
			NDindex[i+I1][t]=counter1;
			counter1+=1;
		}
	}
	
	for (int p=0;p!=P;++p){
		PAN[p]=0;
		PBAset[p].resize(BN);
		for (int b=0;b!=BN;++b){
			PBAN[p][b]=0;
		}
		PINAset[p].resize(ND);
		for (int n=0;n!=ND;++n){
			PINAN[p][n]=0;
		}
		POUTAset[p].resize(ND);
		for (int n=0;n!=ND;++n){
			POUTAN[p][n]=0;
		}
	}
	for (int n=0;n!=ND;++n){
		INAN[n]=0;
		OUTAN[n]=0;
	}
	for (int p1=0;p1!=P;++p1){
		PSAN[p1]=0;  
	}
	counter1=0;
	for (int i1=0;i1!=I1;++i1){
		for (int i2=0;i2!=I2;++i2){
			for (int n=0;n!=RLBN[i1][i2];++n){
				int thebid=LBset[i1][i2][n];
				for (int m=0;m!=SHN[thebid];++m){ 
					int theship=SHset[thebid][m];
					int thestart=SHsts[theship];
					int theend=SHets[theship];
					int nd1=NDindex[i1][thestart];
					int nd2=NDindex[i2+I1][theend];
					int thep=TPH[thestart];  
					int theendp=TPH[theend];  
					PAN[thep]+=1; 
					PAset[thep].resize(PAN[thep],counter1); 
					PBAN[thep][thebid]+=1;  
					PBAset[thep][thebid].resize(PBAN[thep][thebid],counter1); 
					INAN[nd2]+=1;  
					INAset[nd2].resize(INAN[nd2],counter1);
					PINAN[thep][nd2]+=1;
					PINAset[thep][nd2].resize(PINAN[thep][nd2],counter1);  
					for (int cp1=thep; cp1<theendp; ++cp1){
						if (SPmark[cp1][nd2]==-1){
							SPmark[cp1][nd2]=SPANN[cp1];
							SPANN[cp1]+=1;
							SPANset[cp1].resize(SPANN[cp1],nd2); 
						}
					}
					OUTAN[nd1]+=1; 
					OUTAset[nd1].resize(OUTAN[nd1],counter1);
					POUTAN[thep][nd1]+=1;
					POUTAset[thep][nd1].resize(POUTAN[thep][nd1],counter1);		
					if (theendp>thep){
						for (int pp=thep;pp<theendp;++pp){
						     PSAN[pp]+=1; 
							 PSAset[pp].resize(PSAN[pp],counter1); 
						}
					}
					counter1+=1;
					ARCND[0].resize(counter1,nd1); 
					ARCND[1].resize(counter1,nd2); 
					ARCPN.resize(counter1,thep); 
					ARCPS.resize(counter1,PAN[thep]-1); 
					COST3.resize(counter1,c3[thebid]);  
				}
			}
		}
	}
	for (int i1=0;i1!=I1;++i1){
		for (int t1=0;t1!=T;++t1){
			for (int i2=0;i2!=I2;++i2){
				if (t1+len[i1][i2]<T){
					int t2=t1+len[i1][i2];
					int nd1=NDindex[i1][t1];
					int nd2=NDindex[i2+I1][t2];
					int thep=TPH[t1];  
					int theendp=TPH[t2];  
					PAN[thep]+=1; 
					PAset[thep].resize(PAN[thep],counter1); 
					INAN[nd2]+=1;  
					INAset[nd2].resize(INAN[nd2],counter1);
					PINAN[thep][nd2]+=1;
					PINAset[thep][nd2].resize(PINAN[thep][nd2],counter1);
					for (int cp1=thep; cp1<theendp; ++cp1){
						if (SPmark[cp1][nd2]==-1){
							SPmark[cp1][nd2]=SPANN[cp1];
							SPANN[cp1]+=1;
							SPANset[cp1].resize(SPANN[cp1],nd2); 
						}
					}
					OUTAN[nd1]+=1; 
					OUTAset[nd1].resize(OUTAN[nd1],counter1);
					POUTAN[thep][nd1]+=1;
					POUTAset[thep][nd1].resize(POUTAN[thep][nd1],counter1);
					if (theendp>thep){
						for (int pp=thep;pp<theendp;++pp){
							PSAN[pp]+=1; 
							PSAset[pp].resize(PSAN[pp],counter1); 
						}
					}
					counter1+=1;
					ARCND[0].resize(counter1,nd1); 
					ARCND[1].resize(counter1,nd2); 
					ARCPN.resize(counter1,thep); 
					ARCPS.resize(counter1,PAN[thep]-1); 
					COST3.resize(counter1,c4[i1][i2]);  
				}
			}
		}
	}
	ARC=counter1;
	for (int pp=0;pp!=P-1;++pp){  
		for (int nn=0;nn!=SPANN[pp];++nn){
			int thend=SPANset[pp][nn];
			int thesite=Nodeset[0][thend];
			int thetime=Nodeset[1][thend];
			double thedem=0; 
			int startp=pp+1; 
			int endp=TPH[thetime]; 
			for (int crtp=startp;crtp<=endp;++crtp){ 
				 double minitemp=-LL;  
			     int sttime=PTset[crtp][0]; 
			     int edtime=PTset[crtp][PT-1]; 
				 if (edtime>thetime){
					 edtime=thetime;
				 }
				 for (int sn=0;sn!=SN;++sn){
				     double snmini=0;  
					 for (int tt=sttime;tt<=edtime;++tt){
						 int thepseq=tt-sttime;
					     snmini+=sup[thesite][crtp][sn][thepseq];
					 }
					 if (minitemp<snmini){
					 minitemp=snmini;
					 }
				 }
				 thedem+=minitemp;
			}
			MINIdmd[pp].resize(nn+1,thedem);  
		}
	}

}

void Averagesetup(){
	for (int p=0;p!=P;++p){
		for (int i=0;i!=I;++i){
			for (int pt=0;pt!=PTN[p];++pt){
				int thetime=PTset[p][pt];
				int thend=NDindex[i][thetime];
				double thedem=0;
				for (int sn=0;sn!=SN;++sn){
					thedem+=sup[i][p][sn][pt]*PRO[p][sn];
				}
				AVGD[thend]=thedem;
			}
		}
	}
}

void Cplexsetup () { 
	for (int b=0;b!=BN;++b){
		mfb[b]=frt[b]*SHN[b];
	}
	for (int p=0;p!=P;++p){
		for (int i=0;i!=I;++i){
			for (int pt=0;pt!=PTN[p];++pt){
				PNDindex[i][pt]=PNN[p]; 
				int thetime=PTset[p][pt];
				PNN[p]+=1;
				int thend=NDindex[i][thetime];
				PNset[p].resize(PNN[p],thend);
			}
		}
	}
	for (int p=0;p!=P;++p){
		WND[p]=0;
		for (int i=0;i!=I2;++i){
			for (int pt=0;pt!=PTN[p];++pt){
				WNDindex[i][pt]=WND[p];
				WND[p]+=1;
			}
		}
	}
	for (int p=0;p!=P;++p){
		PSAindex[p].resize(ARC,-1);
		for (int nn=0;nn!=PSAN[p];++nn){
		     int thearc=PSAset[p][nn];
			 PSAindex[p][thearc]=nn; 
		}
	}
	Arcflow.resize(ARC,0);
}

void Parasetup(){
	for (int nn=0;nn!=TTRD;++nn){
		NSOL[nn]=625;
	}
	for (int nn=0;nn!=TTRD;++nn){
		RSCE[nn].resize(NSOL[nn]);
		for (int rn=0;rn!=NSOL[nn];++rn){
			RSCE[nn][rn].resize(P);
			for (int pp=0;pp!=P;++pp){
				RSCE[nn][rn][pp]=random(SN);
			}
		}
		VAlsol[nn].resize(NSOL[nn],0);
		for (int nn=0;nn!=TTRD;++nn){
			PD[nn].resize(ND);
		    PInven[nn].resize(P-1);
			PUVval[nn].resize(P-1);
			PArcflow[nn].resize(P-1);
			for (int p=0;p!=P-1;++p){
				PInven[nn][p].resize(I,0);
				PUVval[nn][p].resize(I,0);
				PArcflow[nn][p].resize(PSAN[p],0);
			}
		} 
	}
}

void SOcutsetup(){  
	for (int p1=0;p1!=P;++p1){
		SOARCN[p1]=0;
		SOARCind[p1].resize(ARC,-1);
		SONDN[p1]=0;
		SONDind[p1].resize(ND,-1);
		int aconter=0;
		int nconter=0;
		for (int p2=p1+1;p2!=P;++p2){ 
		     SOARCN[p1]+=PAN[p2]; 
			 for (int nn=0;nn!=PAN[p2];++nn){
			      int thearc=PAset[p2][nn]; 
			      SOARCind[p1][thearc]=aconter;
				  aconter+=1;
				  SOARCset[p1].resize(aconter,thearc);
			 }
			 SONDN[p1]+=PNN[p2]; 
			 for (int nn=0;nn!=PNN[p2];++nn){
			      int thenode=PNset[p2][nn];
				  SONDind[p1][thenode]=nconter;
			      nconter+=1;
				  SONDset[p1].resize(nconter,thenode);
			 }
		}
	}
}

void inputData () { 
	Netsetup();
    Averagesetup();
	Cplexsetup();
	Parasetup();
	SOcutsetup();
}

double CPLEX () {
	double theresult;
	ILOSTLBEGIN
	IloEnv env; 
	IloCplex cplex(env);
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVarArray x(env, BN, 0, 1, ILOINT); 
	IloNumVarArray y(env, BN, 0, IloInfinity, ILOFLOAT); 
	IloNumVarArray z(env, ARC, 0, IloInfinity, ILOFLOAT); 
	IloNumVarArray u(env, ND, 0, IloInfinity, ILOFLOAT); 
	IloNumVarArray v(env, ND, 0, IloInfinity, ILOFLOAT); 
	try{
		
		IloInt b,a,p,n;   
		
		IloModel model(env);
		{   
			model.add(IloMinimize(env, obj)); 
		} 
		{   
			IloExpr Z1(env); 
			IloExpr Z2(env); 
			IloExpr Z3(env); 
			IloExpr Z4(env); 
			for (b=0;b!=BN;++b){
				Z1+=mfb[b]*y[b];  
			}
			for (n=0;n!=ND;++n){
				Z2+=COST1[n]*u[n]; 
			}
			for (n=0;n!=ND;++n){
				Z3+=COST2[n]*v[n]; 
			}
			for (a=0;a!=ARC;++a){
				Z4+=COST3[a]*z[a];  
			}
			model.add(obj==Z1+Z2+Z3+Z4);
			Z1.end();
			Z2.end();
			Z3.end();
			Z4.end();
		}
		
		{   
			for (b=0;b!=BN;++b){
				model.add(y[b]>=lbcap[b]*x[b]);
			}
		}
		{   
			for (b=0;b!=BN;++b){
				model.add(y[b]<=ubcap[b]*x[b]);
			}
		}
		{   
			for (p=0;p!=P;++p){
				for (b=0;b!=BN;++b){
					for (n=0;n!=PBAN[p][b];++n){
						int thearc=PBAset[p][b][n];
						model.add(z[thearc]<=y[b]);
					}
				}
			}
		}
		{   
			for (int t1=1;t1!=T;++t1){
				int t2=t1-1;
				for (int i=0;i!=I1;++i){  
					int nd1=NDindex[i][t1]; 
					int nd2=NDindex[i][t2];
					IloExpr sum1(env);
					for (n=0;n!=OUTAN[nd1];++n){
						int thearc=OUTAset[nd1][n];
						sum1+=z[thearc];
					}
					model.add(u[nd1]==u[nd2]+AVGD[nd1]+v[nd2]-v[nd1]-sum1); 
					sum1.end();
				}
			}
		}
		{   
			for (int i=0;i!=I1;++i){  
				int nd=NDindex[i][0]; 
				IloExpr sum1(env);
				for (n=0;n!=OUTAN[nd];++n){
					int thearc=OUTAset[nd][n];
					sum1+=z[thearc];
				}
				model.add(u[nd]==IQ[nd]+AVGD[nd]-v[nd]-sum1); 
				sum1.end();
			}
		}
		{   
			for (int t1=1;t1!=T;++t1){
				int t2=t1-1;
				for (int i=0;i!=I2;++i){  
					int nd1=NDindex[i+I1][t1]; 
					int nd2=NDindex[i+I1][t2];
					IloExpr sum1(env);
					for (n=0;n!=INAN[nd1];++n){
						int thearc=INAset[nd1][n];
						sum1+=z[thearc];
					}
					model.add(u[nd1]==u[nd2]+AVGD[nd1]-v[nd2]+v[nd1]+sum1); 
					sum1.end();
				}
			}
		}
		{   
			for (int i=0;i!=I2;++i){  
				int nd=NDindex[i+I1][0]; 
				IloExpr sum1(env);
				for (n=0;n!=INAN[nd];++n){
					int thearc=INAset[nd][n];
					sum1+=z[thearc];
				}
				model.add(u[nd]==IQ[nd]+AVGD[nd]+v[nd]+sum1); 
				sum1.end();
			}
		}
		{   
			for (n=0;n!=ND;++n){
				model.add(u[n]<=TQ[n]);
			}
		}
		
		{
			cplex.extract(model);	
			
			cplex.setParam(IloCplex::Threads, 1);  
			cplex.solve();
			cout<<"The model gap is"<<cplex.getMIPRelativeGap()*100<<endl;
			theresult=cplex.getObjValue();
			cout<<"The result is"<<theresult<<endl;			
			cout<<"Bid Selection Results:"<<endl;
			double IHcost=0;  
			double TRcost=0; 
			for (b=0;b!=BN;++b){
				double theres=cplex.getValue(x[b]);
				if (theres>0.99){
					double theamount=cplex.getValue(y[b]);
					Bcap[b]=theamount;  
					cout<<"bid-"<<b<<"is used at"<<theamount<<endl;
					bidcost+=mfb[b]*theamount;
				}
			}
			cout<<"bidcost="<<bidcost<<endl;
		}
		for (a=0;a!=ARC;++a){
			double theflow=cplex.getValue(z[a]);
			if (theflow>0.0001){
				Arcflow[a]=theflow;}
		}
	}
	catch(IloException& ex){
		cerr << ex << endl;
	}
	catch(...){
		cerr << "Error..." << endl;
	}
	env.end();
	return theresult;
}

double TRP02 (int SCE, int PH){   
	double result=0;
	int TASN=0;  
	if (PH>0){
		TASN=PSAN[PH-1]; 
	}
	IloEnv env; 
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVarArray z(env, PAN[PH], 0, IloInfinity, ILOFLOAT); 
	IloNumVarArray u(env, PNN[PH], 0, IloInfinity, ILOFLOAT); 
	IloNumVarArray v(env, PNN[PH], 0, IloInfinity, ILOFLOAT); 
	IloNumVarArray ty(env, BN,-IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray tu(env, I, -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray tv(env, I, -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray tz(env, TASN, -IloInfinity, IloInfinity, ILOFLOAT); 
	try{
		
		IloModel model(env);
		{   
			model.add(IloMinimize(env, obj)); 
		} 
		{   
			IloExpr Z2(env); 
			IloExpr Z3(env); 
			IloExpr Z4(env); 
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				Z2+=COST1[thend]*u[nn];
			}
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				Z3+=COST2[thend]*v[nn];
			}
			for (int nn=0;nn!=PAN[PH];++nn){
				int thearc=PAset[PH][nn];
				Z4+=COST3[thearc]*z[nn];
			}
			model.add(obj==Z2+Z3+Z4);
			Z2.end();
			Z3.end();
			Z4.end();
		}
		
		{   
			for (int b=0;b!=BN;++b){
				for (int n=0;n!=PBAN[PH][b];++n){
					int thearc=PBAset[PH][b][n];
					int theind=ARCPS[thearc]; 
					model.add(z[theind]<=ty[b]);
				}
			}
		}
		{   
			for (int nn=1;nn!=PTN[PH];++nn){ 
				for (int i=0;i!=I1;++i){  
					int ndind1=PNDindex[i][nn];  
					int ndind2=PNDindex[i][nn-1]; 
					int nd1=PNset[PH][ndind1]; 
					
					IloExpr sum1(env);
					for (int n=0;n!=POUTAN[PH][nd1];++n){
						int thearc=POUTAset[PH][nd1][n];
						int theind=ARCPS[thearc]; 
						sum1+=z[theind];
					}
					
					model.add(u[ndind1]+v[ndind1]==u[ndind2]+v[ndind2]+PD[SCE][nd1]-sum1); 
					sum1.end();
				}
			}
		}
		{   
			for (int i=0;i!=I1;++i){
				int ndind1=PNDindex[i][0]; 
				int nd1=PNset[PH][ndind1]; 
				IloExpr sum1(env);
				for (int n=0;n!=POUTAN[PH][nd1];++n){
					int thearc=POUTAset[PH][nd1][n];
					int theind=ARCPS[thearc]; 
					sum1+=z[theind];
				}
				model.add(u[ndind1]+v[ndind1]==tu[i]+tv[i]+PD[SCE][nd1]-sum1);  
				sum1.end();
			}
		}
		{   
			for (int nn=1;nn!=PTN[PH];++nn){ 
				for (int i=0;i!=I2;++i){   
					int ndind1=PNDindex[i+I1][nn]; 
					int ndind2=PNDindex[i+I1][nn-1]; 
					int nd1=PNset[PH][ndind1]; 
					IloExpr tsum1(env); 
					if (PH>0){
						for (int n=0;n!=INAN[nd1];++n){
							int thearc=INAset[nd1][n];
							int theind=PSAindex[PH-1][thearc]; 
							if (theind!=-1){
								tsum1+=tz[theind];}
						}
					}
					IloExpr sum1(env);
					for (int n=0;n!=PINAN[PH][nd1];++n){
						int thearc=PINAset[PH][nd1][n];
						int theind=ARCPS[thearc]; 
						sum1+=z[theind];
					}
					model.add(u[ndind1]-v[ndind1]==u[ndind2]-v[ndind2]+PD[SCE][nd1]+tsum1+sum1); 
					tsum1.end();
					sum1.end();
				}
			}
		}
		{   
			for (int i=0;i!=I2;++i){
				int ndind1=PNDindex[i+I1][0]; 
				int nd1=PNset[PH][ndind1]; 
				IloExpr tsum1(env); 
				if (PH>0){
					for (int n=0;n!=INAN[nd1];++n){
						int thearc=INAset[nd1][n];
						int theind=PSAindex[PH-1][thearc]; 
						if (theind!=-1){
							tsum1+=tz[theind];}
					}
				}
				IloExpr sum1(env);
				for (int n=0;n!=PINAN[PH][nd1];++n){
					int thearc=PINAset[PH][nd1][n];
					int theind=ARCPS[thearc]; 
					sum1+=z[theind];
				}
				model.add(u[ndind1]-v[ndind1]==tu[i+I1]-tv[i+I1]+PD[SCE][nd1]+tsum1+sum1); 
				tsum1.end();
				sum1.end();
			}
		}
		{ 
			if (PH<P-1){  
				for (int nn=0;nn!=SPANN[PH];++nn){
					int thenode=SPANset[PH][nn];      
					int thesite=Nodeset[0][thenode];  
					int lastndind=PNDindex[thesite][PTN[PH]-1]; 
					int thetime=Nodeset[1][thenode];  
					IloExpr sum1(env);
					for (int tt=PTset[PH+1][0]; tt<=thetime;++tt){
						int crtnd=NDindex[thesite][tt]; 
						if (PH>0){
							for (int an=0;an!=PSAN[PH-1];++an){  
								int thearc=PSAset[PH-1][an]; 
								if (ARCND[1][thearc]==crtnd){ 
									sum1+=tz[an];
								}
							}
						}
						for (int an=0;an!=PINAN[PH][crtnd];++an){
							int thearc=PINAset[PH][crtnd][an];
							int theind=ARCPS[thearc]; 
							sum1+=z[theind];
						}
					}
					model.add(u[lastndind]-v[lastndind]+sum1+MINIdmd[PH][nn]<=TQ[thenode]);
					sum1.end();
				}
			}
		}
		{   
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				model.add(u[nn]<=TQ[thend]);
			}
		}
		{   
			for (int b=0;b!=BN;++b){
				model.add(ty[b]==Bcap[b]); 
			}
		}
		{   
			for (int i=0;i!=I;++i){
				if (PH==0){
					model.add(tu[i]==iniv[i]); 
				}
				if (PH>0){
					model.add(tu[i]==PInven[SCE][PH-1][i]); 
				}
			}
		}
		{   
			for (int i=0;i!=I;++i){
				if (PH==0){
					model.add(tv[i]==0); 
				}
				if (PH>0){
					model.add(tv[i]==PUVval[SCE][PH-1][i]); 
				}
			}
		}
		{   
			if (PH>0){
				for (int nn=0;nn!=PSAN[PH-1];++nn){  
					model.add(tz[nn]==PArcflow[SCE][PH-1][nn]);
				}
			}
		}
		
		{
			IloCplex cplex(env);
			cplex.extract(model);
			cplex.setParam(IloCplex::Threads, 1);  
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			cplex.solve();
			double theobj=cplex.getObjValue();
			transcost=theobj; 
			result=transcost; 
			
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				double theval=cplex.getValue(u[nn]);
				STivvol[SCE]+=theval;
				STivcost[SCE]+=theval*COST1[thend]; 
			}
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				double theval=cplex.getValue(v[nn]);
				STbgvol[SCE]+=theval;
				STbgcost[SCE]+=theval*COST2[thend]; 
			}
			for (int nn=0;nn!=PAN[PH];++nn){
				int thearc=PAset[PH][nn];
				double theval=cplex.getValue(z[nn]);
				if (COST3[thearc]<0.001){
					STbidvol[SCE]+=theval;
				}
				if (COST3[thearc]>0.001){
					STspotvol[SCE]+=theval;
					STspotcost[SCE]+=theval*COST3[thearc];
				}
			}
			for (int b=0;b!=BN;++b){
				for (int n=0;n!=PBAN[PH][b];++n){
					int thearc=PBAset[PH][b][n];
					int theind=ARCPS[thearc]; 
					double theval=cplex.getValue(z[theind]);
					STbidvol[SCE]+=theval;
				}
			}
			if (PH<P-1){
				for (int i=0;i!=I;++i){
					int theind=PNDindex[i][PTN[PH]-1]; 
					double theu=cplex.getValue(u[theind]);
					PInven[SCE][PH][i]=theu;
					double thev=cplex.getValue(v[theind]);
					PUVval[SCE][PH][i]=thev;
				}
				for (int nn=0;nn!=PSAN[PH];++nn){
					int thearc=PSAset[PH][nn]; 
					int thep=ARCPN[thearc]; 
					double thef=0;
					if (thep==PH){  
						int theind=ARCPS[thearc];  
						thef=cplex.getValue(z[theind]);
					}
					else {    
						int theind=PSAindex[PH-1][thearc]; 
						thef=cplex.getValue(tz[theind]);
					}
					PArcflow[SCE][PH][nn]=thef; 
				}
			}
		}
	}
	catch(IloException& ex){
		cerr << ex << endl;
	}
	catch(...){
		cerr << "Error..." << endl;
	}
	env.end();
	return result;
}

double SCECAl(int TD, int SIND){
	for (int p=0;p!=P;++p){
        int thess=RSCE[TD][SIND][p];
		for (int n=0;n!=PTN[p];++n){  
			int thetime=PTset[p][n];
			for (int i=0;i!=I;++i){
				int thend=NDindex[i][thetime];
				PD[TD][thend]=sup[i][p][thess][n]; 
			}
		}
	}
	double theupper=bidcost;
	double scetran=0;
	for (int p=0;p!=P-1;++p){ 
		for (int i=0;i!=I;++i){
			PInven[TD][p][i]=0;
			PUVval[TD][p][i]=0;
		}
		PArcflow[TD][p].resize(PSAN[p],0);
	}
	for (int pp=0;pp!=P;++pp){
		double transcost=TRP02(TD,pp);
		scetran+=transcost; 
	}
	theupper+=scetran;
	return theupper;
}

void Upperbound(int M){ 
	cout<<"bidcost="<<bidcost<<endl;
	vector<double> Fres(M);
	UB=0;
#pragma omp parallel sections num_threads(16)
	{
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			
			for (int mm=0;mm!=NSOL[0];++mm){
				cout<<"mm="<<mm*TTRD<<endl;
				double theres=SCECAl(0,mm);
				VAlsol[0][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[1];++mm){
				double theres=SCECAl(1,mm);
				VAlsol[1][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[2];++mm){
				double theres=SCECAl(2,mm);
				VAlsol[2][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[3];++mm){
				double theres=SCECAl(3,mm);
				VAlsol[3][mm]=theres;
			}
		}			
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[4];++mm){
				double theres=SCECAl(4,mm);
				VAlsol[4][mm]=theres;
			}
		}			
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[5];++mm){
				double theres=SCECAl(5,mm);
				VAlsol[5][mm]=theres;
			}
		}	
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[6];++mm){
				double theres=SCECAl(6,mm);
				VAlsol[6][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[7];++mm){
				double theres=SCECAl(7,mm);
				VAlsol[7][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[8];++mm){
				double theres=SCECAl(8,mm);
				VAlsol[8][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[9];++mm){
				double theres=SCECAl(9,mm);
				VAlsol[9][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[10];++mm){
				double theres=SCECAl(10,mm);
				VAlsol[10][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[11];++mm){
				double theres=SCECAl(11,mm);
				VAlsol[11][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[12];++mm){
				double theres=SCECAl(12,mm);
				VAlsol[12][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[13];++mm){
				double theres=SCECAl(13,mm);
				VAlsol[13][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[14];++mm){
				double theres=SCECAl(14,mm);
				VAlsol[14][mm]=theres;
			}
		}
#pragma omp section
		{
			if (MUT<int(omp_get_thread_num())){
				MUT=int(omp_get_thread_num());
			}
			for (int mm=0;mm!=NSOL[15];++mm){
				double theres=SCECAl(15,mm);
				VAlsol[15][mm]=theres;
			}
		}
}	
	int STSN=0;
	for (int nn=0;nn!=TTRD;++nn){
		for (int mm=0;mm!=NSOL[nn];++mm){
			UB+=VAlsol[nn][mm];
			Fres[STSN]=VAlsol[nn][mm];
			STSN+=1;
		}
	}
	cout<<"STSN="<<STSN<<endl;
	double aveub=double(UB)/double(M);
	double totaldev=0;
	for (int nn=0;nn!=M;++nn){
		totaldev+=(double(Fres[nn]-aveub))*(double(Fres[nn]-aveub));
	}
	double avedev=(totaldev)/(M-1);
	cout<<"aveub="<<aveub<<endl;
	cout<<"avedev="<<avedev<<endl;
	UB=aveub+1.96*pow(avedev,0.5)/pow(double(M),0.5);
	stgap=double(UB-LB)/double(LB);
	stgap*=100;
	for (int b=0;b!=BN;++b){
		TTbidcap+=Bcap[b]*SHN[b];
	}
	for (int nn=0;nn!=TTRD;++nn){
			TTspotcost+=STspotcost[nn];
			TTivcost+=STivcost[nn];
			TTbgcost+=STbgcost[nn];
			TTbidvol+=STbidvol[nn];
			TTspotvol+=STspotvol[nn];
			TTivvol+=STivvol[nn];
			TTbgvol+=STbgvol[nn];
	}
	TTspotcost=TTspotcost/double(M);
	TTivcost=TTivcost/double(M);
	TTbgcost=TTbgcost/double(M);
	TTbidvol=TTbidvol/double(M);
	TTspotvol=TTspotvol/double(M);
	TTivvol=TTivvol/double(M);
	TTbgvol=TTbgvol/double(M);
	TTbirate=TTbirate/double(M);
}

void Enumeration(){
	cout<<"bidcost="<<bidcost<<endl;
	double Paradoable[TTRD];
	int search[P];    
	for (int p=0;p!=P-1;++p){
		search[p]=0;   
	}
	double PUB[TTRD];
	for (int tt=0;tt!=TTRD;++tt){ 
		PUB[tt]=0;
	}
	int STSN=0; 
	int doable=1;
	while (doable==1){
		int thes=0;
		while (thes<SN){ 
			for (int tt=0;tt!=TTRD;++tt){ 
				Paradoable[tt]=0;
			}
			for (int tt=0;tt!=TTRD;++tt){  
				Paradoable[tt]=1;
				search[P-1]=thes; 
				for (int p=0;p!=P-1;++p){ 
					for (int n=0;n!=PTN[p];++n){  
						int thetime=PTset[p][n];
						for (int i=0;i!=I;++i){
							int thend=NDindex[i][thetime];
							PD[tt][thend]=sup[i][p][search[p]][n]; 
						}
					}
				}
				for (int n=0;n!=PTN[P-1];++n){  
					int thetime=PTset[P-1][n];
					for (int i=0;i!=I;++i){
						int thend=NDindex[i][thetime];
						PD[tt][thend]=sup[i][P-1][thes][n];
					}
				}
				for (int p=0;p!=P-1;++p){ 
					for (int i=0;i!=I;++i){
						PInven[tt][p][i]=0;
						PUVval[tt][p][i]=0;
					}
					PArcflow[tt][p].resize(PSAN[p],0);
				}
				STSN+=1;
				thes+=1;
				if (thes>=SN){
					break;
				}
			}
			
#pragma omp parallel sections num_threads(16)
			{
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[0]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(0,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[0]+=theupper;
					}
				}
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[1]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(1,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[1]+=theupper;
					}
				}
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[2]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(2,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[2]+=theupper;
					}
				}
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[3]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(3,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[3]+=theupper;
					}
				}			
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[4]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(4,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[4]+=theupper;
					}
				}			
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[5]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(5,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[5]+=theupper;
					}
				}	
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[6]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(6,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[6]+=theupper;
					}
				}	
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[7]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(7,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[7]+=theupper;
					}
				}	

#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[8]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(8,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[8]+=theupper;
					}
				}	

#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[9]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(9,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[9]+=theupper;
					}
				}	

#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[10]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(10,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[10]+=theupper;
					}
				}	

#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[11]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(11,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[11]+=theupper;
					}
				}	
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[12]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(12,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[12]+=theupper;
					}
				}	
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[13]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(13,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[13]+=theupper;
					}
				}	
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[14]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(14,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[14]+=theupper;
					}
				}	
#pragma omp section
				{
					if (MUT<int(omp_get_thread_num())){
						MUT=int(omp_get_thread_num());
					}
					if (Paradoable[15]==1){
						double theupper=bidcost;
						double scetran=0;
						for (int pp=0;pp!=P;++pp){
							double transcost=TRP02(15,pp);
							scetran+=transcost; 
						}
						theupper+=scetran;
						PUB[15]+=theupper;
					}
				}	
			
	}
			cout<<"STSN="<<STSN<<endl;
		}
       
		doable=0;
		for (int p=P-2;p>=0;--p){ 
			if (search[p]<SN-1){
				search[p]+=1;
				doable=1;
				for (int pp=p+1;pp<P-1;++pp){
					search[pp]=0; 
				}
				break;
			}
		}
	}
	UB=0;
	for (int tt=0;tt!=TTRD;++tt){ 
		UB+=PUB[tt];
	}
	double aveub=double(UB)/double(STSN);
	UB=aveub;
	cout<<"aveub para="<<aveub<<endl;
	stgap=double(aveub-LB)/double(LB);
	stgap*=100;
	for (int b=0;b!=BN;++b){
		TTbidcap+=Bcap[b]*SHN[b];
	}
	for (int nn=0;nn!=TTRD;++nn){
		TTspotcost+=STspotcost[nn];
		TTivcost+=STivcost[nn];
		TTbgcost+=STbgcost[nn];
		TTbidvol+=STbidvol[nn];
		TTspotvol+=STspotvol[nn];
		TTivvol+=STivvol[nn];
		TTbgvol+=STbgvol[nn];
	}
	TTspotcost=TTspotcost/double(STSN);
	TTivcost=TTivcost/double(STSN);
	TTbgcost=TTbgcost/double(STSN);
	TTbidvol=TTbidvol/double(STSN);
	TTspotvol=TTspotvol/double(STSN);
	TTivvol=TTivvol/double(STSN);
	TTbgvol=TTbgvol/double(STSN);
	TTbirate=TTbirate/double(STSN);
}

void ORISDDP( ){
    struct timeval p4_start;   
	struct timeval p4_end;   
	gettimeofday(&p4_start,NULL);
	LB=CPLEX ();
	gettimeofday(&p4_end,NULL);
	t_p4=((p4_end.tv_sec - p4_start.tv_sec)*1000000+(p4_end.tv_usec - p4_start.tv_usec))/double(1000000);
	if (P<=7){
		ALLSN=pow(double(SN),double(P));}
	else {
		ALLSN=999999;}
	cout<<"ALLSN="<<ALLSN<<endl;
	if (ALLSN<=ESTV){ 
		cout<<"Do enumeration"<<endl;
		Enumeration();}
	if (ALLSN>ESTV){
		cout<<"Do sampling"<<endl;
		Upperbound(ESTV);
	}
	upperbound=UB;
	totalgap=stgap;
}

void documentation(){
	cout <<"succeeded" << endl;
	cout <<"Threads Used" <<MUT<<endl;
	cout <<"LB=" <<LB<<endl;
	if (ALLSN<=ESTV){ 
		cout<<"Real UB"<<endl;}
	if (ALLSN>ESTV){
		cout<<"Statistical UB"<<endl;
	}
	cout<<"upperbound="<<upperbound<<endl;
	double thegap=(UB-LB)/LB*100;
	cout <<"GAP=" <<thegap<<endl;
	cout <<"Time=" <<t_p4<<endl;
	cout<<"bidcost="<<bidcost<<endl;
	for (int bb=0;bb!=BN;++bb){
		if (Bcap[bb]>1){
			cout<<"bid["<<bb<<"] purchased at"<<Bcap[bb]<<endl;
		}
	}
	cout<<"Tbidcap="<<TTbidcap<<endl;
	cout<<"Tbidcost="<<bidcost<<endl;
	cout<<"Tspotcost="<<TTspotcost<<endl;
	cout<<"Tivcost="<<TTivcost<<endl;
	cout<<"Tbgcost="<<TTbgcost<<endl;
	cout<<"Tbidvol="<<TTbidvol<<endl;
	cout<<"Tspotvol="<<TTspotvol<<endl;
	cout<<"Tivvol="<<TTivvol<<endl;
	cout<<"Tbgvol="<<TTbgvol<<endl;
}

int main (int argc, char **argv) {
	srand(2);
	srand((int)time(NULL));
	inputData();
	ORISDDP();
    documentation();
	return 0;
}
