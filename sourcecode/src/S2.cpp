#include <ilcplex/ilocplex.h>   /* cplex */
#include <stdio.h>              /* printf, scanf, puts, NULL */
#include <stdlib.h>             /* srand, rand */
#include <sys/time.h>              /* time */             
#include <fstream>  
#include <string>
#include <math.h>
#include <sstream>
#include <algorithm>            /* random_shuffle */
#include <vector>
#include "Avgminmax02.h"
#include "Seqinsertion.h"
#define random(x) (rand()%x)
namespace patch {template < typename T > std::string to_string( const T& n ){std::ostringstream stm; stm << n; return stm.str() ;}}

using namespace std;

//input data starts here
int Dev = 10;// Deviation
extern const int CNN = 1;// Case Index
extern const int P = 3;// Number of phases
extern const int PT = 6;// Number of periods in a phase
extern const int T = 18;// Number of periods
extern const int I1 =1;// Number of suppliers
extern const int I2 = 8;// Number of customers
extern const int I = 9;// Number of sites
extern const int BN = 72;//Number of bids
extern const int LBN = 9;//Number of bids on each arc
extern const int SPN = 339;//Number of shipments
extern const int SN = 10;//Number of scenarios at each stage
extern const int MBSN = 9;//Number of shipments in a bid
// tree parameters:
double PRO[P][SN]={{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1},{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1},{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1}}; // the probobility of each scenario at each phase 
// timing parameters: 
int PTN[P]={6,6,6}; // number of periods in each phase 
int PTset[P][PT]={{0,1,2,3,4,5},{6,7,8,9,10,11},{12,13,14,15,16,17}}; // set of periods in each phase 
// distribution parameters: 
int sup[I][P][SN][PT]={{{{11408,0,0,0,0,0},{11615,0,0,0,0,0},{11400,0,0,0,0,0},{11421,0,0,0,0,0},{11558,0,0,0,0,0},{11534,0,0,0,0,0},{11594,0,0,0,0,0},{11595,0,0,0,0,0},{11726,0,0,0,0,0},{11553,0,0,0,0,0}},{{11316,0,0,0,0,0},{11677,0,0,0,0,0},{11516,0,0,0,0,0},{11630,0,0,0,0,0},{11392,0,0,0,0,0},{11482,0,0,0,0,0},{11441,0,0,0,0,0},{11563,0,0,0,0,0},{11761,0,0,0,0,0},{11582,0,0,0,0,0}},{{11667,0,0,0,0,0},{11389,0,0,0,0,0},{11415,0,0,0,0,0},{11419,0,0,0,0,0},{11460,0,0,0,0,0},{11422,0,0,0,0,0},{11466,0,0,0,0,0},{11689,0,0,0,0,0},{11472,0,0,0,0,0},{11566,0,0,0,0,0}}},{{{-185,-221,-193,-187,-199,-209},{-195,-197,-201,-185,-207,-209},{-189,-183,-185,-209,-219,-201},{-199,-183,-185,-203,-209,-219},{-195,-195,-223,-193,-217,-201},{-221,-195,-207,-197,-213,-201},{-203,-217,-215,-183,-215,-187},{-217,-211,-197,-217,-213,-187},{-213,-223,-189,-207,-221,-217},{-215,-215,-185,-201,-183,-189}},{{-213,-191,-199,-183,-185,-215},{-195,-219,-213,-203,-209,-215},{-195,-215,-209,-219,-221,-221},{-205,-189,-223,-189,-189,-205},{-203,-199,-201,-199,-187,-215},{-203,-205,-201,-199,-215,-213},{-203,-191,-191,-205,-207,-183},{-193,-187,-223,-213,-197,-223},{-223,-199,-223,-197,-191,-223},{-191,-203,-195,-193,-219,-217}},{{-215,-223,-217,-221,-195,-205},{-203,-215,-199,-193,-187,-213},{-209,-207,-221,-207,-183,-189},{-215,-221,-189,-217,-191,-207},{-195,-189,-197,-223,-199,-205},{-191,-207,-203,-207,-211,-187},{-223,-211,-215,-193,-183,-195},{-211,-201,-213,-221,-207,-217},{-207,-201,-223,-197,-219,-191},{-189,-205,-209,-205,-223,-195}}},{{{-172,-156,-148,-169,-163,-150},{-161,-167,-171,-159,-167,-167},{-147,-147,-174,-169,-155,-156},{-159,-161,-155,-163,-159,-150},{-161,-169,-148,-177,-153,-175},{-145,-172,-177,-156,-156,-163},{-153,-177,-145,-166,-177,-169},{-150,-164,-159,-150,-156,-172},{-158,-172,-167,-159,-164,-167},{-147,-150,-159,-177,-171,-172}},{{-150,-164,-167,-155,-171,-161},{-156,-175,-163,-163,-177,-145},{-148,-155,-167,-166,-155,-155},{-175,-145,-167,-175,-145,-150},{-153,-171,-169,-159,-161,-172},{-172,-147,-156,-158,-156,-172},{-147,-150,-148,-147,-164,-163},{-174,-167,-147,-163,-148,-164},{-158,-166,-158,-167,-175,-145},{-172,-169,-150,-156,-172,-174}},{{-164,-158,-148,-156,-171,-158},{-150,-172,-156,-172,-145,-172},{-151,-164,-163,-175,-169,-161},{-161,-161,-158,-151,-163,-163},{-174,-159,-161,-148,-164,-156},{-167,-153,-159,-163,-145,-147},{-155,-169,-158,-153,-175,-163},{-158,-169,-151,-175,-171,-163},{-175,-156,-150,-174,-167,-167},{-166,-163,-155,-159,-153,-150}}},{{{-231,-229,-243,-208,-238,-236},{-240,-210,-254,-249,-213,-247},{-226,-245,-222,-215,-222,-238},{-238,-249,-238,-224,-213,-245},{-213,-240,-245,-217,-231,-224},{-210,-231,-213,-226,-219,-243},{-210,-245,-231,-236,-245,-236},{-231,-224,-217,-236,-243,-240},{-240,-217,-229,-222,-247,-252},{-247,-254,-231,-213,-238,-238}},{{-222,-217,-247,-229,-238,-215},{-213,-217,-219,-238,-231,-245},{-229,-247,-231,-215,-210,-252},{-233,-219,-236,-222,-233,-254},{-215,-254,-222,-233,-224,-222},{-224,-252,-238,-217,-215,-231},{-233,-226,-210,-243,-233,-219},{-224,-254,-249,-233,-226,-217},{-245,-240,-249,-231,-247,-254},{-249,-240,-233,-222,-236,-226}},{{-245,-247,-247,-226,-245,-231},{-236,-236,-249,-245,-249,-213},{-226,-252,-252,-231,-247,-219},{-224,-233,-249,-226,-217,-252},{-233,-213,-240,-245,-229,-215},{-224,-247,-222,-217,-233,-236},{-236,-245,-236,-229,-208,-217},{-226,-215,-229,-245,-252,-247},{-229,-243,-215,-219,-245,-217},{-238,-252,-233,-229,-224,-219}}},{{{-226,-217,-217,-202,-234,-221},{-234,-215,-200,-221,-228,-221},{-215,-228,-204,-221,-217,-239},{-195,-202,-239,-228,-224,-224},{-226,-217,-197,-204,-206,-239},{-219,-200,-213,-234,-195,-224},{-237,-206,-234,-197,-204,-204},{-226,-210,-239,-226,-195,-213},{-202,-224,-226,-210,-210,-215},{-226,-208,-239,-217,-208,-230}},{{-202,-197,-232,-226,-208,-197},{-215,-210,-213,-221,-221,-210},{-197,-221,-230,-200,-224,-221},{-200,-221,-232,-195,-208,-232},{-208,-234,-195,-208,-195,-197},{-239,-206,-221,-206,-224,-226},{-213,-228,-217,-230,-206,-228},{-234,-206,-213,-239,-206,-224},{-230,-208,-195,-234,-200,-210},{-215,-221,-221,-206,-213,-200}},{{-202,-204,-224,-232,-197,-195},{-228,-215,-219,-206,-239,-217},{-224,-224,-232,-195,-217,-237},{-237,-215,-195,-204,-232,-202},{-219,-239,-237,-208,-195,-200},{-230,-221,-221,-221,-221,-206},{-210,-202,-200,-195,-210,-221},{-237,-234,-195,-237,-215,-221},{-197,-224,-195,-197,-200,-237},{-234,-197,-213,-219,-224,-221}}},{{{-279,-293,-247,-271,-242,-261},{-287,-263,-247,-293,-255,-277},{-261,-245,-287,-277,-253,-279},{-250,-239,-261,-261,-274,-277},{-269,-285,-255,-242,-245,-269},{-253,-279,-279,-282,-287,-285},{-285,-293,-279,-245,-261,-290},{-282,-287,-282,-266,-247,-287},{-271,-293,-277,-274,-269,-269},{-250,-263,-282,-250,-253,-263}},{{-247,-239,-258,-255,-239,-271},{-271,-293,-282,-239,-263,-258},{-293,-239,-258,-242,-271,-290},{-287,-293,-282,-279,-245,-290},{-250,-242,-261,-285,-282,-269},{-263,-261,-274,-245,-271,-255},{-263,-239,-282,-255,-277,-261},{-253,-274,-255,-258,-271,-239},{-271,-261,-290,-285,-269,-245},{-274,-250,-247,-247,-279,-261}},{{-287,-269,-290,-266,-285,-245},{-245,-293,-239,-242,-266,-261},{-255,-258,-285,-271,-239,-239},{-266,-279,-250,-250,-282,-255},{-250,-239,-293,-274,-263,-263},{-274,-247,-253,-242,-282,-266},{-239,-282,-266,-245,-287,-261},{-255,-287,-290,-282,-263,-250},{-271,-253,-293,-253,-266,-255},{-261,-287,-242,-253,-274,-258}}},{{{-300,-260,-272,-252,-263,-280},{-255,-274,-294,-283,-274,-258},{-252,-291,-260,-266,-294,-277},{-252,-255,-300,-288,-269,-286},{-297,-308,-286,-283,-302,-288},{-258,-274,-269,-297,-283,-288},{-300,-302,-269,-260,-266,-252},{-291,-263,-308,-308,-277,-286},{-272,-300,-277,-258,-294,-294},{-305,-258,-302,-288,-305,-283}},{{-294,-308,-297,-277,-291,-274},{-277,-258,-302,-308,-300,-291},{-294,-302,-252,-283,-263,-266},{-305,-272,-272,-305,-269,-294},{-255,-263,-255,-297,-272,-297},{-283,-305,-263,-274,-260,-252},{-297,-258,-283,-286,-274,-269},{-272,-291,-266,-305,-277,-260},{-308,-274,-288,-297,-308,-269},{-266,-277,-263,-297,-300,-283}},{{-291,-297,-288,-252,-258,-294},{-288,-258,-252,-255,-288,-294},{-252,-252,-252,-263,-291,-302},{-260,-302,-258,-255,-300,-280},{-255,-302,-266,-269,-286,-305},{-300,-286,-300,-269,-302,-300},{-274,-308,-294,-291,-294,-266},{-260,-305,-286,-286,-252,-283},{-263,-272,-288,-302,-294,-274},{-258,-255,-291,-272,-286,-308}}},{{{-300,-260,-294,-266,-291,-308},{-288,-286,-283,-302,-277,-252},{-272,-258,-252,-294,-300,-288},{-260,-305,-277,-260,-269,-294},{-308,-294,-286,-288,-266,-297},{-260,-286,-280,-274,-260,-272},{-308,-272,-291,-286,-277,-272},{-266,-291,-272,-288,-252,-269},{-274,-294,-255,-288,-302,-286},{-252,-283,-260,-286,-302,-308}},{{-283,-255,-263,-260,-255,-255},{-291,-283,-305,-269,-269,-305},{-283,-288,-308,-266,-252,-288},{-300,-300,-255,-269,-300,-266},{-288,-280,-277,-283,-288,-272},{-280,-302,-269,-269,-286,-297},{-277,-283,-302,-297,-308,-255},{-252,-291,-308,-305,-269,-297},{-294,-277,-288,-288,-255,-291},{-294,-286,-272,-277,-297,-297}},{{-302,-272,-255,-305,-297,-255},{-280,-291,-269,-274,-286,-258},{-291,-274,-280,-272,-252,-269},{-283,-255,-272,-260,-300,-255},{-300,-266,-280,-286,-266,-302},{-294,-272,-258,-260,-288,-277},{-263,-258,-277,-255,-266,-280},{-291,-280,-297,-277,-266,-269},{-260,-283,-288,-269,-305,-258},{-252,-269,-305,-302,-294,-302}}},{{{-291,-266,-277,-263,-258,-260},{-280,-308,-302,-274,-302,-283},{-260,-291,-255,-280,-305,-277},{-297,-288,-269,-263,-286,-277},{-272,-283,-260,-277,-272,-260},{-291,-305,-288,-291,-291,-272},{-272,-255,-302,-280,-283,-302},{-277,-274,-305,-283,-258,-283},{-300,-305,-283,-302,-286,-252},{-269,-255,-283,-252,-297,-291}},{{-260,-294,-283,-283,-286,-305},{-291,-286,-302,-252,-302,-294},{-297,-263,-308,-266,-255,-286},{-274,-297,-260,-305,-300,-269},{-294,-305,-263,-272,-283,-263},{-274,-263,-297,-297,-266,-280},{-274,-305,-266,-291,-302,-252},{-294,-274,-288,-260,-283,-297},{-305,-266,-280,-274,-308,-302},{-305,-266,-286,-300,-263,-302}},{{-288,-266,-294,-283,-294,-308},{-300,-274,-258,-255,-274,-260},{-263,-283,-297,-263,-280,-277},{-308,-291,-263,-272,-260,-280},{-260,-308,-263,-280,-255,-286},{-263,-302,-258,-252,-280,-258},{-263,-308,-294,-277,-308,-308},{-286,-308,-291,-263,-255,-297},{-302,-255,-286,-286,-263,-288},{-294,-300,-263,-308,-260,-297}}}}; // supply (positive for supply and negative for demand) under different scenarios 
int iniv[I]={2070,261,23,150,170,304,400,300,160}; // initial inventory 
int ubiv[I]={2740,522,414,300,341,608,800,600,320}; // the upper bounds of inventory 
int len[I1][I2]={{3,6,2,2,3,3,2,1}}; // the traveling time between two sites 
// cost parameters: 
double c1[I]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // inventory holding cost 
double c2[I]={0.05,7.19,13.32,4.32,2.64,5.16,6.85,4.63,1.82}; // cost of leftover and unmet 
double c3[BN]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; // transportation cost for capacity contracts 
double c4[I1][I2]={{6.54,12.11,3.93,2.4,4.69,6.23,4.21,1.65}}; // cost of spot market 
// bid parameters: 
int RLBN[I1][I2]={{9,9,9,9,9,9,9,9}}; // number of bids on each arc 
int LBset[I1][I2][LBN]={{{0,1,2,3,4,5,6,7,8},{9,10,11,12,13,14,15,16,17},{18,19,20,21,22,23,24,25,26},{27,28,29,30,31,32,33,34,35},{36,37,38,39,40,41,42,43,44},{45,46,47,48,49,50,51,52,53},{54,55,56,57,58,59,60,61,62},{63,64,65,66,67,68,69,70,71}}}; // set of bids on each arc 
double frt[BN]={5.23,5.23,5.23,4.58,4.58,4.58,3.92,3.92,3.92,9.69,9.69,9.69,8.48,8.48,8.48,7.27,7.27,7.27,3.14,3.14,3.14,2.75,2.75,2.75,2.36,2.36,2.36,1.92,1.92,1.92,1.68,1.68,1.68,1.44,1.44,1.44,3.75,3.75,3.75,3.28,3.28,3.28,2.81,2.81,2.81,4.98,4.98,4.98,4.36,4.36,4.36,3.74,3.74,3.74,3.37,3.37,3.37,2.95,2.95,2.95,2.53,2.53,2.53,1.32,1.32,1.32,1.16,1.16,1.16,0.99,0.99,0.99}; // freight rate of each bid 
int lbcap[BN] = { 75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226}; // capacity lower bound of each bid  
int ubcap[BN]={150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300}; // capacity upper bound of each bid 
// shipment parameters: 
int SHN[BN]={8,4,3,7,4,3,8,4,3,6,3,2,5,3,2,6,3,2,7,4,3,7,4,3,8,4,3,7,4,3,7,4,3,8,4,3,8,4,3,7,4,3,8,4,3,8,4,3,7,4,3,7,4,3,7,4,3,7,4,3,8,4,3,9,4,3,9,4,3,8,4,3}; // number of shipment for each bid 
int SHset[BN][MBSN]={{0,1,2,3,4,5,6,7},{8,9,10,11},{12,13,14},{15,16,17,18,19,20,21},{22,23,24,25},{26,27,28},{29,30,31,32,33,34,35,36},{37,38,39,40},{41,42,43},{44,45,46,47,48,49},{50,51,52},{53,54},{55,56,57,58,59},{60,61,62},{63,64},{65,66,67,68,69,70},{71,72,73},{74,75},{76,77,78,79,80,81,82},{83,84,85,86},{87,88,89},{90,91,92,93,94,95,96},{97,98,99,100},{101,102,103},{104,105,106,107,108,109,110,111},{112,113,114,115},{116,117,118},{119,120,121,122,123,124,125},{126,127,128,129},{130,131,132},{133,134,135,136,137,138,139},{140,141,142,143},{144,145,146},{147,148,149,150,151,152,153,154},{155,156,157,158},{159,160,161},{162,163,164,165,166,167,168,169},{170,171,172,173},{174,175,176},{177,178,179,180,181,182,183},{184,185,186,187},{188,189,190},{191,192,193,194,195,196,197,198},{199,200,201,202},{203,204,205},{206,207,208,209,210,211,212,213},{214,215,216,217},{218,219,220},{221,222,223,224,225,226,227},{228,229,230,231},{232,233,234},{235,236,237,238,239,240,241},{242,243,244,245},{246,247,248},{249,250,251,252,253,254,255},{256,257,258,259},{260,261,262},{263,264,265,266,267,268,269},{270,271,272,273},{274,275,276},{277,278,279,280,281,282,283,284},{285,286,287,288},{289,290,291},{292,293,294,295,296,297,298,299,300},{301,302,303,304},{305,306,307},{308,309,310,311,312,313,314,315,316},{317,318,319,320},{321,322,323},{324,325,326,327,328,329,330,331},{332,333,334,335},{336,337,338}}; // set of shipments in each bid 
int SHsts[SPN]={0,2,4,6,8,10,12,14,2,6,10,14,1,7,13,2,4,6,8,10,12,14,2,6,10,14,0,6,12,0,2,4,6,8,10,12,14,2,6,10,14,0,6,12,0,2,4,6,8,10,2,6,10,1,7,2,4,6,8,10,0,4,8,2,8,1,3,5,7,9,11,1,5,9,1,7,2,4,6,8,10,12,14,0,4,8,12,1,7,13,2,4,6,8,10,12,14,1,5,9,13,0,6,12,0,2,4,6,8,10,12,14,2,6,10,14,0,6,12,2,4,6,8,10,12,14,1,5,9,13,1,7,13,2,4,6,8,10,12,14,0,4,8,12,1,7,13,0,2,4,6,8,10,12,14,2,6,10,14,0,6,12,0,2,4,6,8,10,12,14,1,5,9,13,0,6,12,2,4,6,8,10,12,14,2,6,10,14,2,8,14,0,2,4,6,8,10,12,14,0,4,8,12,2,8,14,0,2,4,6,8,10,12,14,0,4,8,12,1,7,13,2,4,6,8,10,12,14,2,6,10,14,1,7,13,2,4,6,8,10,12,14,0,4,8,12,2,8,14,2,4,6,8,10,12,14,0,4,8,12,1,7,13,2,4,6,8,10,12,14,1,5,9,13,0,6,12,1,3,5,7,9,11,13,15,1,5,9,13,1,7,13,0,2,4,6,8,10,12,14,16,1,5,9,13,0,6,12,0,2,4,6,8,10,12,14,16,1,5,9,13,1,7,13,2,4,6,8,10,12,14,16,1,5,9,13,2,8,14}; // set of start times of each shipment 
int SHets[SPN]={3,5,7,9,11,13,15,17,5,9,13,17,4,10,16,5,7,9,11,13,15,17,5,9,13,17,3,9,15,3,5,7,9,11,13,15,17,5,9,13,17,3,9,15,6,8,10,12,14,16,8,12,16,7,13,8,10,12,14,16,6,10,14,8,14,7,9,11,13,15,17,7,11,15,7,13,4,6,8,10,12,14,16,2,6,10,14,3,9,15,4,6,8,10,12,14,16,3,7,11,15,2,8,14,2,4,6,8,10,12,14,16,4,8,12,16,2,8,14,4,6,8,10,12,14,16,3,7,11,15,3,9,15,4,6,8,10,12,14,16,2,6,10,14,3,9,15,2,4,6,8,10,12,14,16,4,8,12,16,2,8,14,3,5,7,9,11,13,15,17,4,8,12,16,3,9,15,5,7,9,11,13,15,17,5,9,13,17,5,11,17,3,5,7,9,11,13,15,17,3,7,11,15,5,11,17,3,5,7,9,11,13,15,17,3,7,11,15,4,10,16,5,7,9,11,13,15,17,5,9,13,17,4,10,16,5,7,9,11,13,15,17,3,7,11,15,5,11,17,4,6,8,10,12,14,16,2,6,10,14,3,9,15,4,6,8,10,12,14,16,3,7,11,15,2,8,14,3,5,7,9,11,13,15,17,3,7,11,15,3,9,15,1,3,5,7,9,11,13,15,17,2,6,10,14,1,7,13,1,3,5,7,9,11,13,15,17,2,6,10,14,2,8,14,3,5,7,9,11,13,15,17,2,6,10,14,3,9,15}; // set of end times of each shipment 

//input data ends here

int LL=999999;
double sddperr=-100;  // set it to a negative value means that we do not allow stablize to terminate
double uptime_p1=1560*P;
double uptime_p2=1200*P;
double pdlgap=1; 
double pdlgap02=1;
double timlim01=600;
double timlim02=1200;
int SAMN=16;
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
double UB=LL;
double LB=-LL;
double bidcost;
double thebidlb; 
double transcost;
double stgap=100; 

double Bcap[BN];    
vector<vector<vector<double>>> Inven;   
vector<vector<vector<double>>> UVval;   
vector<vector<vector<double>>> Arcflow; 

int SPANN[P];                    
vector<vector<int>> SPANset(P-1);  
vector<vector<double>> MINIdmd(P-1);   
int SPmark[P][ND];     

int C4N[P];  
int C5N[P];  
int C6N=I1;  
int C7N[P];  
int C8N=I2;  
int C9N[P];  
int C10N[P];  
int C28N=BN; 
int C29N=I;  
int C30N=I;  
int C31N[P];  

vector<vector<int>> Ac4ind(P); 
vector<vector<int>> Ac5ind(P);
vector<vector<int>> Ac6ind(P); 
vector<vector<int>> Ac7ind(P); 
vector<vector<int>> TAc7ind(P); 
vector<vector<int>> Ac8ind(P); 
vector<vector<int>> TAc8ind(P); 
vector<vector<int>> TAc31ind(P); 
vector<vector<vector<int>>> Arcc10ind(P-1); 
vector<vector<vector<int>>> TArcc10ind(P-1); 
vector<vector<int>> Nc5ind1(P);
vector<vector<int>> Nc5ind2(P);
vector<vector<int>> Nc6ind(P);
vector<vector<int>> TNc6ind(P);
vector<vector<int>> Nc7ind1(P);
vector<vector<int>> Nc7ind2(P);
vector<vector<int>> Nc8ind(P);
vector<vector<int>> TNc8ind(P);
vector<vector<int>> Nc9ind(P);
vector<vector<int>> TNc29ind(P);
vector<vector<int>> TNc30ind(P);
vector<vector<vector<int>>> Ndc10ind(P-1);  

int BidCN=0;  
vector<double> Bidcutval; 
vector<vector<double>> Bidcutset(BN); 

int TRCN[P]; 
vector<vector<double>> TRcutval(P-1); 
vector<vector<vector<double>>> TRzcutset(P-1); 
vector<vector<vector<double>>> TRucutset(P-1); 
vector<vector<vector<double>>> TRvcutset(P-1); 
vector<vector<vector<double>>> TRtycutset(P-1); 

int TSN=16;
vector<double> trpr;  

double SInven[P][I];   
double SUVval[P][I];   
vector<vector<double>> SArcflow(P); 


double upperbound;
double lb_total;
double lb_p1;
double lb_p2;
double lb_p3;

double totalgap;
double t_total; 
double t_p1;
double t_p2;
double t_p3;
double t_p4;

int itrn_total=0;
int itrn_p1=0;
int itrn_p2=0;
int tbcn_total=0;
int tbcn_p1=0;
int tbcn_p2=0;
int ALLSN; 
int RTSN;  

double PD[ND];  //supply and demand
double PInven[P-1][I];   // the inventory of each node at the end of each phase
double PUVval[P-1][I];   // the overflow and the unmet at the end of each phase
vector<vector<double>> PArcflow(P-1); // the flow on each arc that could affect the next stages 

double PTRtycutset[BN];
double PTRucutset[I];
double PTRvcutset[I];
vector<vector<double>> PTRzcutset(P);
double PTRcutval;
double PBidcutset[BN];  // the coefficient for each bid in the benders cut
double PBidcutval;


int Xval[BN];    
double Bgap;
double Bidf;      
vector<vector<double>> Ffval(P);     

int STTSN[P];  
vector<vector<double>> STtrpr(P);  
vector<vector<int>> STGSSN(P); 
vector<vector<vector<int>>> STGPSN(P); 
vector<vector<vector<double>>> STGSpro(P); 
vector<vector<vector<vector<vector<double>>>>> STGsup(P); 
vector<vector<vector<int>>> STGSSCE(P);  
vector<vector<int>> STGPATHN(P); 
vector<vector<vector<int>>> STGSFSN(P); 
vector<vector<vector<vector<int>>>> STGSFSset(P); 
vector<vector<vector<double>>> STD(P);  

int SOARCN[P]; 
int SONDN[P];  
vector<vector<int>> SOARCset(P); 
vector<vector<int>> SONDset(P);   
vector<vector<int>> SOARCind (P); 
vector<vector<int>> SONDind (P); 
int MUT=0; 

double TTspotcost=0;  
double TTivcost=0;    
double TTbgcost=0;    
double TTbidcap=0;   
double TTbidvol=0;    
double TTspotvol=0;   
double TTivvol=0;    
double TTbgvol=0;    
double TTivvol_s=0;    
double TTbgvol_s=0;    
double TTivvol_d=0;    
double TTbgvol_d=0;    
double TTbirate=0;   


string num2str(double i)
{    stringstream ss;
ss<<i;
return ss.str();}

ofstream mycout01("//scratch//lxwu//2024//HSDDP//Results//D"+num2str(Dev)+"_P"+num2str(P)+"_SN"+num2str(SN)+"_CN"+num2str(CNN)+".txt");

void Subtreesetup(){ 
	for (int thep=0;thep!=P;++thep){
		STGSSN[thep].resize(P,1);
		int ADL=2;
		if (thep>0){
			ADL=0;
		}
		int endp=thep+ADL;
		if (endp>P-1){
			endp=P-1;
		}
		int thetsn=1;
		int ASN=2;
		for (int pp=thep;pp<=endp;++pp){
			STGSSN[thep][pp]=ASN;  
			thetsn=thetsn*ASN;
		}
		STTSN[thep]=thetsn;
		STGPSN[thep].resize(P);
		for (int nn=thep;nn!=P;++nn){
			
			STGPSN[thep][nn].resize(STGSSN[thep][nn],0);
			int TBSN=floor(double(SN)/double(STGSSN[thep][nn]));
			for (int gg=0;gg!=STGSSN[thep][nn];++gg){
				STGPSN[thep][nn][gg]=TBSN;
			}
			int Toadd=SN-TBSN*STGSSN[thep][nn];
			for (int aa=0;aa!=Toadd;++aa){
				STGPSN[thep][nn][aa]+=1;
			}
		}
		STGsup[thep].resize(I);
		for (int i=0;i!=I;++i){
			STGsup[thep][i].resize(P);
			for (int p=thep;p!=P;++p){
				STGsup[thep][i][p].resize(STGSSN[thep][p]);
				for (int sn=0;sn!=STGSSN[thep][p];++sn){
					STGsup[thep][i][p][sn].resize(PT,0);
				}
			}
		}
		STGSpro[thep].resize(P);
		for (int p=thep;p!=P;++p){
			STGSpro[thep][p].resize(STGSSN[thep][p],0);
		}
		STGPATHN[thep].resize(P);
		STGPATHN[thep][thep]=STGSSN[thep][thep]; 
		STGSFSN[thep].resize(P);
		STGSFSN[thep][thep].resize(STGPATHN[thep][thep],0);
		STGSFSset[thep].resize(P);
		STGSFSset[thep][thep].resize(STGPATHN[thep][thep]);
		for (int p=thep+1;p!=P;++p){
			STGPATHN[thep][p]=STGPATHN[thep][p-1]*STGSSN[thep][p];  
			STGSFSN[thep][p].resize(STGPATHN[thep][p],0); 
			STGSFSset[thep][p].resize(STGPATHN[thep][p]);  
		}
		STtrpr[thep].resize(STTSN[thep],0); 
		int search[P];
		for (int p=thep;p!=P;++p){
			search[p]=0;
		}
		int doable=1;
		int counter=0;
		STGSSCE[thep].resize(P);
		while (doable==1){
			for (int s=0;s!=STGSSN[thep][P-1];++s){  
				search[P-1]=s;  
				STGSSCE[thep][P-1].resize(counter+1,s); 
				for (int p=thep;p!=P-1;++p){ 
					STGSSCE[thep][p].resize(counter+1,search[p]); 
				}
				int thenode=search[thep];
				STGSFSN[thep][thep][thenode]+=1;
				STGSFSset[thep][thep][thenode].resize(STGSFSN[thep][thep][thenode],counter); 
				for (int p=thep+1;p!=P;++p){
					thenode=thenode*STGSSN[thep][p]+search[p]; 
					STGSFSN[thep][p][thenode]+=1;
					STGSFSset[thep][p][thenode].resize(STGSFSN[thep][p][thenode],counter); 
				}
				counter+=1; 
			}
			
			doable=0;
			for (int p=P-2;p>=thep;--p){ 
				if (search[p]<STGSSN[thep][p]-1){
					search[p]+=1;
					doable=1;
					for (int pp=p+1;pp<P;++pp){
						search[pp]=0; 
					}
					break;
				}
			}
		}
	}
}

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
}

void SDDPsetup(){ 
	for (int bb=0;bb!=BN;++bb){
		Bcap[bb]=0;
	}
	for (int p=0;p!=P-1;++p){
		PArcflow[p].resize(PSAN[p],0);
	}
	for(int pp=0;pp!=P;++pp){
		PTRzcutset[pp].resize(PSAN[pp],0); // the coefficient for z variables
	}
	for (int p=0;p!=P-1;++p){
	     TRCN[p]=0;
		 TRucutset[p].resize(I);       
		 TRvcutset[p].resize(I);       
		 TRtycutset[p].resize(BN);     
		 TRzcutset[p].resize(PSAN[p]); 
	}
	TRCN[P-1]=0; 
	D.resize(TSN);
	for (int nn=0;nn!=TSN;++nn){
		D[nn].resize(ND);} 
	Inven.resize(TSN);
	UVval.resize(TSN);
	Arcflow.resize(TSN);
	for (int sn=0;sn!=TSN;++sn){
		Inven[sn].resize(P-1);
		UVval[sn].resize(P-1);
		Arcflow[sn].resize(P-1);
		for (int p=0;p!=P-1;++p){
			Inven[sn][p].resize(I,0);
			UVval[sn][p].resize(I,0);
			Arcflow[sn][p].resize(PSAN[p],0);
		}
	}
}

void Dualsetup(){ 
	for (int p=0;p!=P-1;++p){ 
		Arcc10ind[p].resize(SPANN[p]);
		Ndc10ind[p].resize(SPANN[p]);
		TArcc10ind[p].resize(SPANN[p]);
		int TSAN=0;
		if (p>0){
			TSAN=PSAN[p-1];
		}
		for (int nn=0;nn!=SPANN[p];++nn){
			Arcc10ind[p][nn].resize(PAN[p],-1);
			TArcc10ind[p][nn].resize(TSAN,-1);
			Ndc10ind[p][nn].resize(PNN[p],-1);
		}
	}
	for (int p=0;p!=P;++p){  
	     Ac4ind[p].resize(PAN[p],-1);
	     Ac5ind[p].resize(PAN[p],-1);
		 Ac6ind[p].resize(PAN[p],-1);
		 Ac7ind[p].resize(PAN[p],-1);
		 int TSAN=0;
		 if (p>0){
			 TSAN=PSAN[p-1];
		 }
		 TAc7ind[p].resize(TSAN,-1);
		 Ac8ind[p].resize(PAN[p],-1);
		 TAc8ind[p].resize(TSAN,-1);
		 TAc31ind[p].resize(TSAN,-1); 
	}
	for (int p=0;p!=P;++p){  
		Nc5ind1[p].resize(PNN[p],-1);
		Nc5ind2[p].resize(PNN[p],-1);
		Nc6ind[p].resize(PNN[p],-1);
		TNc6ind[p].resize(I,-1);
		Nc7ind1[p].resize(PNN[p],-1);
		Nc7ind2[p].resize(PNN[p],-1);
		Nc8ind[p].resize(PNN[p],-1);
		TNc8ind[p].resize(I,-1);
		Nc9ind[p].resize(PNN[p],-1);
		TNc29ind[p].resize(I,-1);
		TNc30ind[p].resize(I,-1);
	}
	for (int p=0;p!=P;++p){
		C4N[p]=0; 
		for (int b=0;b!=BN;++b){
			for (int n=0;n!=PBAN[p][b];++n){
				int thearc=PBAset[p][b][n];
				int theind=ARCPS[thearc];
				Ac4ind[p][theind]=C4N[p];   
				C4N[p]+=1;
			}
		}
	}
	for (int p=0;p!=P;++p) {   
		C5N[p]=0;
		for (int nn=1;nn!=PTN[p];++nn){ 
			for (int i=0;i!=I1;++i){  
				int ndind1=PNDindex[i][nn];  
				Nc5ind1[p][ndind1]=C5N[p];
				int ndind2=PNDindex[i][nn-1];  
				Nc5ind2[p][ndind2]=C5N[p];  
				int nd1=PNset[p][ndind1]; 
				for (int n=0;n!=POUTAN[p][nd1];++n){
					int thearc=POUTAset[p][nd1][n];
					int theind=ARCPS[thearc]; 
					Ac5ind[p][theind]=C5N[p];
				}
				C5N[p]+=1;
			}
		}
	}
	for (int p=0;p!=P;++p) {  
		for (int i=0;i!=I1;++i){
			int ndind1=PNDindex[i][0]; 
			Nc6ind[p][ndind1]=i;
			TNc6ind[p][i]=i;            
			int nd1=PNset[p][ndind1]; 
			for (int n=0;n!=POUTAN[p][nd1];++n){
				 int thearc=POUTAset[p][nd1][n];
				 int theind=ARCPS[thearc]; 
                 Ac6ind[p][theind]=i; 
			}
		}
	}
	for (int p=0;p!=P;++p) {   
		C7N[p]=0;
		for (int nn=1;nn!=PTN[p];++nn){ 
			for (int i=0;i!=I2;++i){   
				int ndind1=PNDindex[i+I1][nn]; 
				Nc7ind1[p][ndind1]=C7N[p];
				int ndind2=PNDindex[i+I1][nn-1]; 
				Nc7ind2[p][ndind2]=C7N[p];
				int nd1=PNset[p][ndind1]; 
				for (int n=0;n!=PINAN[p][nd1];++n){ 
					int thearc=PINAset[p][nd1][n];
					int theind=ARCPS[thearc]; 
                    Ac7ind[p][theind]=C7N[p];
				}
				if (p>0){
					for (int n=0;n!=INAN[nd1];++n){
						int thearc=INAset[nd1][n];
						
						int theind=PSAindex[p-1][thearc]; 
						
						if (theind!=-1){
							TAc7ind[p][theind]=C7N[p];} 
					}
				}
				C7N[p]+=1;
			}
		}
	}
	for (int p=0;p!=P;++p) {   
		for (int i=0;i!=I2;++i){
			int ndind1=PNDindex[i+I1][0]; 
			Nc8ind[p][ndind1]=i;
			TNc8ind[p][i+I1]=i;
			int nd1=PNset[p][ndind1]; 
			for (int n=0;n!=PINAN[p][nd1];++n){
				int thearc=PINAset[p][nd1][n];
				int theind=ARCPS[thearc]; 
				Ac8ind[p][theind]=i;
			}
			if (p>0){
				for (int n=0;n!=INAN[nd1];++n){
					int thearc=INAset[nd1][n];
					int theind=PSAindex[p-1][thearc]; 
					if (theind!=-1){
						TAc8ind[p][theind]=i;} 
				}
			}
		}
	}
	for (int p=0;p!=P;++p) {   
		C9N[p]=0;
		for (int nn=0;nn!=PNN[p];++nn){
			Nc9ind[p][nn]=C9N[p];
			C9N[p]+=1;
		}
	}
	for (int p=0;p!=P;++p) {   
		for (int i=0;i!=I;++i){
		     TNc29ind[p][i]=i;
		}
	}
	for (int p=0;p!=P;++p) {   
		for (int i=0;i!=I;++i){
			TNc30ind[p][i]=i;
		}
	}
	C31N[0]=0;
	for (int p=0;p!=P;++p) {   
		if (p>0){
			C31N[p]=PSAN[p-1];
			for (int nn=0;nn!=PSAN[p-1];++nn){
				TAc31ind[p][nn]=nn; 
			}
		}
	}
	{ 
		for (int p=0;p!=P-1;++p){  
			C10N[p]=0; 
			for (int nn=0;nn!=SPANN[p];++nn){
				int thenode=SPANset[p][nn];      
				int thesite=Nodeset[0][thenode];  
				int lastndind=PNDindex[thesite][PTN[p]-1]; 
				Ndc10ind[p][nn][lastndind]=C10N[p];  
				int thetime=Nodeset[1][thenode];  
				for (int tt=PTset[p+1][0]; tt<=thetime;++tt){
					int crtnd=NDindex[thesite][tt]; 
					for (int an=0;an!=PINAN[p][crtnd];++an){
						int thearc=PINAset[p][crtnd][an];
						int theind=ARCPS[thearc]; 
						Arcc10ind[p][nn][theind]=C10N[p];  
					}
					if (p>0){
						for (int an=0;an!=PSAN[p-1];++an){  
							int thearc=PSAset[p-1][an]; 
							if (ARCND[1][thearc]==crtnd){ 
								TArcc10ind[p][nn][an]=C10N[p];  
							}
						}
					}
				}
				C10N[p]+=1;
			}
		}
	}
	C10N[P-1]=0; 

}

void Subsampling(){ 
	for (int thep=0;thep!=P;++thep){
		for (int nn=0;nn!=STTSN[thep];++nn){
			STtrpr[thep][nn]=0;
		}
		STD[thep].resize(STTSN[thep]);
		for (int nn=0;nn!=STTSN[thep];++nn){
			STD[thep][nn].resize(ND,0); 
		}
		for (int i=0;i!=I;++i){
			for (int p=thep;p!=P;++p){
				for (int sn=0;sn!=STGSSN[thep][p];++sn){
					for (int t=0;t!=PT;++t){
						STGsup[thep][i][p][sn][t]=0;}
				}
			}
		}
		for (int p=thep;p!=P;++p){
			for (int sn=0;sn!=STGSSN[thep][p];++sn){
				STGSpro[thep][p][sn]=0;}
		}
		
		for (int p=thep;p!=P;++p){
			int Selind[SN];
			int Seq[SN];
			for (int nn=0;nn!=SN;++nn){
				int totalsup=0;
				for (int t=0;t!=PT;++t){
					totalsup+=sup[0][p][nn][t];
				}
				Selind[nn]=totalsup;
				Seq[nn]=nn;
			}
			InsertionSort(Selind,Seq,SN);
			int counter=0;
			for (int gg=0;gg!=STGSSN[thep][p];++gg){  
				int thegn=STGPSN[thep][p][gg];
				double trate=0;
				for (int gsn=counter;gsn!=counter+thegn;++gsn){
					int thesce=Seq[gsn];  
					trate+=PRO[p][thesce]; 
					for (int i=0;i!=I;++i){
						for (int t=0;t!=PT;++t){
							STGsup[thep][i][p][gg][t]+=double(double(sup[i][p][thesce][t])*double(PRO[p][thesce]));
						}
					}
				}
				counter+=thegn;
				STGSpro[thep][p][gg]=trate;
				for (int i=0;i!=I;++i){
					for (int t=0;t!=PT;++t){
						STGsup[thep][i][p][gg][t]=double(double(STGsup[thep][i][p][gg][t])/double(trate));
					}
				}
			}
		}
		
		for (int sn=0;sn!=STTSN[thep];++sn){
			double thepro=1;
			for (int p=thep;p!=P;++p){
				int theps=STGSSCE[thep][p][sn];
				thepro=thepro*STGSpro[thep][p][theps];
				for (int n=0;n!=PTN[p];++n){  
					int thetime=PTset[p][n];
					for (int i=0;i!=I;++i){
						int thend=NDindex[i][thetime];
						STD[thep][sn][thend]=STGsup[thep][i][p][theps][n]; 
					}
				}
			}
			STtrpr[thep][sn]=thepro;
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
	Subtreesetup();
	Netsetup();
	Cplexsetup ();
	Dualsetup();
	SDDPsetup();
	Subsampling();
	SOcutsetup();
}

void Sampling(int M){ 
	vector<double>().swap(trpr); 
	for (int nn=0;nn!=TSN;++nn){
	     vector<double>().swap(D[nn]); 
	}
	double temppro=double(1)/double(M);
	D.resize(M);
	for (int nn=0;nn!=M;++nn){
		trpr.resize(nn+1,temppro);
		D[nn].resize(ND);}
	int STSN=0; 
	while (STSN<M){ 
		for (int p=0;p!=P;++p){
			int thess=random(SN);  
			for (int n=0;n!=PTN[p];++n){  
				int thetime=PTset[p][n];
				for (int i=0;i!=I;++i){
					int thend=NDindex[i][thetime];
					D[STSN][thend]=sup[i][p][thess][n]; 
				}
			}
		}
		STSN+=1; 
	}
	for (int p=0;p!=P;++p){
		Ffval[p].resize(M,0);
	}
}

double BIDP (int record, int LCI) { 
	double result=0;
	IloEnv env; 
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVarArray x(env, BN, 0, 1, ILOFLOAT); 
	IloNumVarArray y(env, BN, 0, IloInfinity, ILOFLOAT); 
	IloNumVar F(env, 0, IloInfinity, ILOFLOAT); 
	IloArray<IloNumVarArray> z(env, STTSN[0]);
	IloArray<IloNumVarArray> u(env, STTSN[0]);
	IloArray<IloNumVarArray> v(env, STTSN[0]);
	for (int sn=0;sn!=STTSN[0];++sn){ 
		z[sn]=IloNumVarArray (env, ARC, 0, IloInfinity, ILOFLOAT); 
		u[sn]=IloNumVarArray (env, ND, 0, IloInfinity, ILOFLOAT); 
		v[sn]=IloNumVarArray (env, ND, 0, IloInfinity, ILOFLOAT); 
	}
	try{
		
		IloModel model(env);
		{   
			model.add(IloMinimize(env, obj)); 
		} 
		{   
			IloExpr Z1(env); 
			for (int b=0;b!=BN;++b){
				Z1+=mfb[b]*y[b];  
			}
			model.add(obj==Z1+F);
			Z1.end();
		}
		
		{ 
			if (LCI==2){ 
				for (int b=0;b!=BN;++b){
					model.add(x[b]==Xval[b]);
				}
			}
		}
		{   
			for (int b=0;b!=BN;++b){
				model.add(y[b]>=lbcap[b]*x[b]);
			}
		}
		{   
			for (int b=0;b!=BN;++b){
				model.add(y[b]<=ubcap[b]*x[b]);
			}
		}
		{    
			for (int n=0;n!=BidCN;++n){
				IloExpr sum1(env);
				for (int bb=0;bb!=BN;++bb){
					sum1+=y[bb]*Bidcutset[bb][n];
				}
				model.add(F>=sum1+Bidcutval[n]);
				sum1.end();
			}
		}
		{    
			IloExpr sum1(env);
			for (int s=0;s!=STTSN[0];++s){
				for (int p=0;p!=P;++p){
					for (int n=0;n!=PNN[p];++n){
						int thend=PNset[p][n]; 
						sum1+=STtrpr[0][s]*COST1[thend]*u[s][thend]; 
						sum1+=STtrpr[0][s]*COST2[thend]*v[s][thend]; 
					}
				}
			}
			for (int s=0;s!=STTSN[0];++s){
				for (int p=0;p!=P;++p){
					for (int n=0;n!=PAN[p];++n){
						int tharc=PAset[p][n]; 
						sum1+=STtrpr[0][s]*COST3[tharc]*z[s][tharc]; 
					}
				}
			}
			model.add(F>=sum1);
			sum1.end();
		}
		{   
			for (int s=0;s!=STTSN[0];++s){
				for (int b=0;b!=BN;++b){
					for (int p=0;p!=P;++p){
						for (int n=0;n!=PBAN[p][b];++n){
							int thearc=PBAset[p][b][n];
							model.add(z[s][thearc]<=y[b]);
						}
					}
				}
			}
		}
		{   
			for (int s=0;s!=STTSN[0];++s){
				for (int t1=1;t1!=T;++t1){
					int t2=t1-1;
					for (int i=0;i!=I1;++i){  
						int nd1=NDindex[i][t1]; 
						int nd2=NDindex[i][t2];
						IloExpr sum1(env);
						for (int n=0;n!=OUTAN[nd1];++n){
							int thearc=OUTAset[nd1][n];
							sum1+=z[s][thearc];
						}
						model.add(u[s][nd1]==u[s][nd2]+STD[0][s][nd1]+v[s][nd2]-v[s][nd1]-sum1); 
						sum1.end();
					}
				}
			}
		}
		{   
			for (int s=0;s!=STTSN[0];++s){
				for (int i=0;i!=I1;++i){  
					int nd=NDindex[i][0]; 
					IloExpr sum1(env);
					for (int n=0;n!=OUTAN[nd];++n){
						int thearc=OUTAset[nd][n];
						sum1+=z[s][thearc];
					}
					model.add(u[s][nd]==IQ[nd]+STD[0][s][nd]-v[s][nd]-sum1); 
					sum1.end();
				}
			}
		}
		{   
			for (int s=0;s!=STTSN[0];++s){
				for (int t1=1;t1!=T;++t1){
					int t2=t1-1;
					for (int i=0;i!=I2;++i){  
						int nd1=NDindex[i+I1][t1]; 
						int nd2=NDindex[i+I1][t2];
						IloExpr sum1(env);
						for (int n=0;n!=INAN[nd1];++n){
							int thearc=INAset[nd1][n];
							sum1+=z[s][thearc];
						}
						model.add(u[s][nd1]==u[s][nd2]+STD[0][s][nd1]-v[s][nd2]+v[s][nd1]+sum1); 
						sum1.end();
					}
				}
			}
		}
		{   
			for (int s=0;s!=STTSN[0];++s){
				for (int i=0;i!=I2;++i){  
					int nd=NDindex[i+I1][0]; 
					IloExpr sum1(env);
					for (int n=0;n!=INAN[nd];++n){
						int thearc=INAset[nd][n];
						sum1+=z[s][thearc];
					}
					model.add(u[s][nd]==IQ[nd]+STD[0][s][nd]+v[s][nd]+sum1); 
					sum1.end();
				}
			}
		}
		{   
			for (int s=0;s!=STTSN[0];++s){
				for (int n=0;n!=ND;++n){
					model.add(u[s][n]<=TQ[n]);
				}
			}
		}
		{   
			for (int p=0;p!=P;++p){
				for (int pn=0;pn!= STGPATHN[0][p];++pn){ 
					for (int sn=0;sn<STGSFSN[0][p][pn]-1;++sn){
						int s1=STGSFSset[0][p][pn][sn];
						for (int mm=sn+1;mm<STGSFSN[0][p][pn];++mm){
							int s2=STGSFSset[0][p][pn][mm];
							for (int an=0;an!=PAN[p];++an){  
								int thearc=PAset[p][an];
								model.add(z[s1][thearc]==z[s2][thearc]);
							}
						}
					}
				}
			}
		}
		{ 
			for (int s=0;s!=STTSN[0];++s){
				for (int PH=0;PH<P-1;++PH){  
					int lasttime=PTset[PH][PTN[PH]-1];
					for (int nn=0;nn!=SPANN[PH];++nn){
						
						int thenode=SPANset[PH][nn];      
						int thesite=Nodeset[0][thenode];  
						int thetime=Nodeset[1][thenode];  
						int lastnd=NDindex[thesite][lasttime];
						IloExpr sum1(env);
						for (int tt=PTset[PH+1][0]; tt<=thetime;++tt){
							int crtnd=NDindex[thesite][tt]; 
							for (int pp=0;pp<=PH;++pp){
								for (int an=0;an!=PINAN[pp][crtnd];++an){
									int thearc=PINAset[pp][crtnd][an];
									sum1+=z[s][thearc];
								}
							}
						}
						model.add(u[s][lastnd]-v[s][lastnd]+sum1+MINIdmd[PH][nn]<=TQ[thenode]);
						sum1.end();
					}
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
			double fval=cplex.getValue(F); 
			Bidf=fval;
			bidcost=theobj-fval;
			result=theobj;
			if (LCI!=2){ 
				thebidlb=theobj;}
			if (record==1){
				for (int b=0;b!=BN;++b){
					Bcap[b]=0; 
					double theamount=cplex.getValue(y[b]);
					Bcap[b]=theamount; 
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

double FLOATBIDP(int record, double timlim) {
	double result = 0;
	IloEnv env;
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVarArray x(env, BN, 0, 1, ILOFLOAT);
	IloNumVarArray y(env, BN, 0, IloInfinity, ILOFLOAT);
	IloNumVar F(env, 0, IloInfinity, ILOFLOAT);
	IloArray<IloNumVarArray> z(env, STTSN[0]);
	IloArray<IloNumVarArray> u(env, STTSN[0]);
	IloArray<IloNumVarArray> v(env, STTSN[0]);
	for (int sn = 0; sn != STTSN[0]; ++sn) {
		z[sn] = IloNumVarArray(env, ARC, 0, IloInfinity, ILOFLOAT);
		u[sn] = IloNumVarArray(env, ND, 0, IloInfinity, ILOFLOAT);
		v[sn] = IloNumVarArray(env, ND, 0, IloInfinity, ILOFLOAT);
	}
	try {

		IloModel model(env);
		{
			model.add(IloMinimize(env, obj));
		}
		{
			IloExpr Z1(env);
			for (int b = 0; b != BN; ++b) {
				Z1 += mfb[b] * y[b];
			}
			model.add(obj == Z1 + F);
			Z1.end();
		}

		{
			for (int b = 0; b != BN; ++b) {
				model.add(y[b] >= lbcap[b] * x[b]);
			}
		}

		{
			for (int b = 0; b != BN; ++b) {
				model.add(y[b] <= ubcap[b] * x[b]);
			}
		}
		{
			for (int n = 0; n != BidCN; ++n) {
				IloExpr sum1(env);
				for (int bb = 0; bb != BN; ++bb) {
					sum1 += y[bb] * Bidcutset[bb][n];
				}
				model.add(F >= sum1 + Bidcutval[n]);
				sum1.end();
			}
		}
		{
			IloExpr sum1(env);
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int p = 0; p != P; ++p) {
					for (int n = 0; n != PNN[p]; ++n) {
						int thend = PNset[p][n];
						sum1 += STtrpr[0][s] * COST1[thend] * u[s][thend];
						sum1 += STtrpr[0][s] * COST2[thend] * v[s][thend];
					}
				}
			}
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int p = 0; p != P; ++p) {
					for (int n = 0; n != PAN[p]; ++n) {
						int tharc = PAset[p][n];
						sum1 += STtrpr[0][s] * COST3[tharc] * z[s][tharc];
					}
				}
			}
			model.add(F >= sum1);
			sum1.end();
		}
		{
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int b = 0; b != BN; ++b) {
					for (int p = 0; p != P; ++p) {
						for (int n = 0; n != PBAN[p][b]; ++n) {
							int thearc = PBAset[p][b][n];
							model.add(z[s][thearc] <= y[b]);
						}
					}
				}
			}
		}
		{
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int t1 = 1; t1 != T; ++t1) {
					int t2 = t1 - 1;
					for (int i = 0; i != I1; ++i) {
						int nd1 = NDindex[i][t1];
						int nd2 = NDindex[i][t2];
						IloExpr sum1(env);
						for (int n = 0; n != OUTAN[nd1]; ++n) {
							int thearc = OUTAset[nd1][n];
							sum1 += z[s][thearc];
						}
						model.add(u[s][nd1] == u[s][nd2] + STD[0][s][nd1] + v[s][nd2] - v[s][nd1] - sum1);
						sum1.end();
					}
				}
			}
		}
		{
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int i = 0; i != I1; ++i) {
					int nd = NDindex[i][0];
					IloExpr sum1(env);
					for (int n = 0; n != OUTAN[nd]; ++n) {
						int thearc = OUTAset[nd][n];
						sum1 += z[s][thearc];
					}
					model.add(u[s][nd] == IQ[nd] + STD[0][s][nd] - v[s][nd] - sum1);
					sum1.end();
				}
			}
		}
		{
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int t1 = 1; t1 != T; ++t1) {
					int t2 = t1 - 1;
					for (int i = 0; i != I2; ++i) {
						int nd1 = NDindex[i + I1][t1];
						int nd2 = NDindex[i + I1][t2];
						IloExpr sum1(env);
						for (int n = 0; n != INAN[nd1]; ++n) {
							int thearc = INAset[nd1][n];
							sum1 += z[s][thearc];
						}
						model.add(u[s][nd1] == u[s][nd2] + STD[0][s][nd1] - v[s][nd2] + v[s][nd1] + sum1);
						sum1.end();
					}
				}
			}
		}
		{
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int i = 0; i != I2; ++i) {
					int nd = NDindex[i + I1][0];
					IloExpr sum1(env);
					for (int n = 0; n != INAN[nd]; ++n) {
						int thearc = INAset[nd][n];
						sum1 += z[s][thearc];
					}
					model.add(u[s][nd] == IQ[nd] + STD[0][s][nd] + v[s][nd] + sum1);
					sum1.end();
				}
			}
		}
		{
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int n = 0; n != ND; ++n) {
					model.add(u[s][n] <= TQ[n]);
				}
			}
		}
		{
			for (int p = 0; p != P; ++p) {
				for (int pn = 0; pn != STGPATHN[0][p]; ++pn) {
					for (int sn = 0; sn < STGSFSN[0][p][pn] - 1; ++sn) {
						int s1 = STGSFSset[0][p][pn][sn];
						for (int mm = sn + 1; mm < STGSFSN[0][p][pn]; ++mm) {
							int s2 = STGSFSset[0][p][pn][mm];
							for (int an = 0; an != PAN[p]; ++an) {
								int thearc = PAset[p][an];
								model.add(z[s1][thearc] == z[s2][thearc]);
							}
						}
					}
				}
			}
		}
		{
			for (int s = 0; s != STTSN[0]; ++s) {
				for (int PH = 0; PH < P - 1; ++PH) {
					int lasttime = PTset[PH][PTN[PH] - 1];
					for (int nn = 0; nn != SPANN[PH]; ++nn) {
						int thenode = SPANset[PH][nn];
						int thesite = Nodeset[0][thenode];
						int thetime = Nodeset[1][thenode];
						int lastnd = NDindex[thesite][lasttime];
						IloExpr sum1(env);
						for (int tt = PTset[PH + 1][0]; tt <= thetime; ++tt) {
							int crtnd = NDindex[thesite][tt];
							for (int pp = 0; pp <= PH; ++pp) {
								for (int an = 0; an != PINAN[pp][crtnd]; ++an) {
									int thearc = PINAset[pp][crtnd][an];
									sum1 += z[s][thearc];
								}
							}
						}
						model.add(u[s][lastnd] - v[s][lastnd] + sum1 + MINIdmd[PH][nn] <= TQ[thenode]);
						sum1.end();
					}
				}
			}
		}

		{
			IloCplex cplex(env);
			cplex.setParam(IloCplex::Param::TimeLimit, timlim);
			cplex.extract(model);
			cplex.setParam(IloCplex::Threads, 1);
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			clock_t t_start01, t_end01;
			t_start01 = clock();
			cplex.solve();
			t_end01 = clock();
			double theobj = cplex.getObjValue();
			double fval = cplex.getValue(F);
			Bidf = fval;
			bidcost = theobj - fval;
			result = theobj;
			if (double(t_end01 - t_start01) / 1000 >= timlim - 1) {
				cout << "time reached" << endl;
				double thegap = cplex.getMIPRelativeGap();
				double thelim = 0.001;
				if (thegap * 100 > thelim) {
					result = result / (1 + thegap);
				}
			}
			if (record == 1) {
				thebidlb = theobj;
				for (int b = 0; b != BN; ++b) {
					int thex = 0;
					Bcap[b] = 0;
					double temp = cplex.getValue(x[b]);
					if (temp >= 0.5) {
						thex = 1;
					}
					Xval[b] = thex;
					if (Xval[b] >= 0.9) {
						double theamount = cplex.getValue(y[b]);
						if (theamount < lbcap[b]) {
							Bcap[b] = lbcap[b];
						}
						if (theamount > ubcap[b]) {
							Bcap[b] = ubcap[b];
						}
						if (theamount >= lbcap[b] && theamount <= ubcap[b]) {
							Bcap[b] = theamount;
						}
						cout << "bid[" << b << "] is selected with amount [" << Bcap[b] << "]" << endl;
					}
				}
			}
		}
	}
	catch (IloException& ex) {
		cerr << ex << endl;
	}
	catch (...) {
		cerr << "Error..." << endl;
	}
	env.end();
	return result;
}

double TRP (int SCE, int PH){   
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
	IloNumVar F(env, 0, IloInfinity, ILOFLOAT); 
	IloArray<IloNumVarArray> sz(env, STTSN[PH+1]);
	IloArray<IloNumVarArray> su(env, STTSN[PH+1]);
	IloArray<IloNumVarArray> sv(env, STTSN[PH+1]);
	for (int sn=0;sn!=STTSN[PH+1];++sn){ 
		sz[sn]=IloNumVarArray (env, SOARCN[PH], 0, IloInfinity, ILOFLOAT); 
		su[sn]=IloNumVarArray (env, SONDN[PH], 0, IloInfinity, ILOFLOAT); 
		sv[sn]=IloNumVarArray (env, SONDN[PH], 0, IloInfinity, ILOFLOAT); 
	}
	try{
		
		IloModel model(env);
		{   
			model.add(IloMinimize(env, obj)); 
		} 
		{   
			IloExpr Z2(env); 
			IloExpr Z3(env); 
			IloExpr Z4(env); 
			IloExpr Z5(env); 
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
			model.add(obj==Z2+Z3+Z4+F);
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
					model.add(u[ndind1]+v[ndind1]==u[ndind2]+v[ndind2]+D[SCE][nd1]-sum1); 
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
				
				model.add(u[ndind1]+v[ndind1]==tu[i]+tv[i]+D[SCE][nd1]-sum1);  
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
					model.add(u[ndind1]-v[ndind1]==u[ndind2]-v[ndind2]+D[SCE][nd1]+tsum1+sum1); 
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
				
				model.add(u[ndind1]-v[ndind1]==tu[i+I1]-tv[i+I1]+D[SCE][nd1]+tsum1+sum1); 
				tsum1.end();
				sum1.end();
			}
		}
		{   
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				model.add(u[nn]<=TQ[thend]);
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
					
					model.add(tu[i]==Inven[SCE][PH-1][i]); 
				}
			}
		}
		{   
			for (int i=0;i!=I;++i){
				if (PH==0){
					model.add(tv[i]==0); 
				}
				if (PH>0){
					model.add(tv[i]==UVval[SCE][PH-1][i]); 
				}
			}
		}
		{   
			if (PH>0){
				for (int nn=0;nn!=PSAN[PH-1];++nn){  
					model.add(tz[nn]==Arcflow[SCE][PH-1][nn]);
				}
			}
		}
		{   
			if (PH<P-1){ 
				
				for (int n=0;n!=TRCN[PH];++n){ 
					IloExpr sum1(env);
					for (int i=0;i!=I;++i){
						double theval=TRucutset[PH][i][n]; 
						int un=PNDindex[i][PTN[PH]-1];     
						sum1+=theval*u[un];
					}
					for (int i=0;i!=I;++i){
						double theval=TRvcutset[PH][i][n]; 
						int un=PNDindex[i][PTN[PH]-1];     
						sum1+=theval*v[un];
					}
					for (int zn=0;zn!=PSAN[PH];++zn){
						double theval=TRzcutset[PH][zn][n]; 
						int thearc=PSAset[PH][zn];   
						int thep=ARCPN[thearc]; 
						if (thep==PH){  
							int theind=ARCPS[thearc];  
							sum1+=theval*z[theind];
						}
						else {    
							int theind=PSAindex[PH-1][thearc];
							sum1+=theval*tz[theind]; 
						}
					}	
					for (int yn=0;yn!=BN;++yn){
						double theval=TRtycutset[PH][yn][n]; 
						sum1+=theval*ty[yn];
					}
					model.add(F>=sum1+TRcutval[PH][n]);
					sum1.end();
				}
			}
		}
		
		{   
			for (int sn=0;sn!=STTSN[PH+1];++sn){ 
				for (int b=0;b!=BN;++b){
					for (int pp=PH+1;pp!=P;++pp){  
						for (int n=0;n!=PBAN[pp][b];++n){
							int thearc=PBAset[pp][b][n];
							int theind=SOARCind[PH][thearc]; 
							model.add(sz[sn][theind]<=ty[b]);
						}
					}
				}
			}
		}
		{   
			for (int sn=0;sn!=STTSN[PH+1];++sn){
				for (int t1=PTset[PH+1][1];t1!=T;++t1){  
					int t2=t1-1; 
					for (int i=0;i!=I1;++i){  
						int nd1=NDindex[i][t1]; 
						int nd2=NDindex[i][t2];
						int ndind1=SONDind[PH][nd1];  
						int ndind2=SONDind[PH][nd2];  
						IloExpr sum1(env);
						for (int n=0;n!=OUTAN[nd1];++n){
							int thearc=OUTAset[nd1][n];
							int theind=SOARCind[PH][thearc]; 
							sum1+=sz[sn][theind];
						}
						model.add(su[sn][ndind1]+sv[sn][ndind1]==su[sn][ndind2]+sv[sn][ndind2]+STD[PH+1][sn][nd1]-sum1); 
						sum1.end();
					}
				}
			}			
		}
		{   
			for (int sn=0;sn!=STTSN[PH+1];++sn){
				int t1=PTset[PH+1][0];  
				for (int i=0;i!=I1;++i){
					int oldndind=PNDindex[i][PTN[PH]-1]; 
					int nd1=NDindex[i][t1]; 
					int ndind=SONDind[PH][nd1]; 
					IloExpr sum1(env);
					for (int n=0;n!=OUTAN[nd1];++n){
						int thearc=OUTAset[nd1][n];
						int theind=SOARCind[PH][thearc]; 
						sum1+=sz[sn][theind];
					}
					
					model.add(su[sn][ndind]+sv[sn][ndind]==u[oldndind]+v[oldndind]+STD[PH+1][sn][nd1]-sum1);  
					sum1.end();
				}
			}
		}
		{   
			for (int sn=0;sn!=STTSN[PH+1];++sn){
				for (int t1=PTset[PH+1][1];t1!=T;++t1){ 
					int t2=t1-1;
					for (int i=0;i!=I2;++i){   
						int nd1=NDindex[i+I1][t1]; 
						int nd2=NDindex[i+I1][t2];
						IloExpr sum1(env);  
						for (int n=0;n!=INAN[nd1];++n){
							int thearc=INAset[nd1][n];
							int thephase=ARCPN[thearc]; 
							if (thephase<PH){ 
								int theind=PSAindex[PH-1][thearc]; 
								sum1+=tz[theind];
							}
							if (thephase==PH){ 
								int theind=ARCPS[thearc]; 
								sum1+=z[theind];
							}
							if (thephase>PH){ 
								int arcind=SOARCind[PH][thearc]; 
								sum1+=sz[sn][arcind];
							}
						}
						int ndind1=SONDind[PH][nd1];  
						int ndind2=SONDind[PH][nd2];  
						model.add(su[sn][ndind1]-sv[sn][ndind1]==su[sn][ndind2]-sv[sn][ndind2]+STD[PH+1][sn][nd1]+sum1); 
						sum1.end();
					}
				}
			}
		}
		{   
			for (int sn=0;sn!=STTSN[PH+1];++sn){
				int t1=PTset[PH+1][0];  
				for (int i=0;i!=I2;++i){   
					int oldndind=PNDindex[i+I1][PTN[PH]-1]; 
					int nd1=NDindex[i+I1][t1]; 
					IloExpr sum1(env);  
					for (int n=0;n!=INAN[nd1];++n){
						int thearc=INAset[nd1][n];
						int thephase=ARCPN[thearc]; 
						if (thephase<PH){ 
							int theind=PSAindex[PH-1][thearc]; 
							sum1+=tz[theind];
						}
						if (thephase==PH){ 
							int theind=ARCPS[thearc]; 
							sum1+=z[theind];
						}
						if (thephase>PH){ 
							int arcind=SOARCind[PH][thearc]; 
							sum1+=sz[sn][arcind];
						}
					}
					int ndind1=SONDind[PH][nd1];  
					model.add(su[sn][ndind1]-sv[sn][ndind1]==u[oldndind]-v[oldndind]+STD[PH+1][sn][nd1]+sum1); 
					sum1.end();
				}
			}
		}
		{   
			for (int sn=0;sn!=STTSN[PH+1];++sn){
				for (int nn=0;nn!=SONDN[PH];++nn){
					int thend=SONDset[PH][nn]; 
					model.add(su[sn][nn]<=TQ[thend]);
				}
			}
		}
		{ 
			for (int sn=0;sn!=STTSN[PH+1];++sn){
				for (int pp=PH+1;pp!=P-1;++pp){
					for (int nn=0;nn!=SPANN[pp];++nn){
						
						int thenode=SPANset[pp][nn];      
						int thesite=Nodeset[0][thenode];  
						int lastnd=NDindex[thesite][PTN[pp]-1]; 
						int thetime=Nodeset[1][thenode];  
						IloExpr sum1(env);
						for (int tt=PTset[pp+1][0]; tt<=thetime;++tt){
							int crtnd=NDindex[thesite][tt]; 
							for (int n=0;n!=INAN[crtnd];++n){
								int thearc=INAset[crtnd][n];
								int thephase=ARCPN[thearc]; 
								if (thephase<PH){ 
									int theind=PSAindex[PH-1][thearc]; 
									sum1+=tz[theind];
								}
								if (thephase==PH){ 
									int theind=ARCPS[thearc]; 
									sum1+=z[theind];
								}
								if (thephase>PH){ 
									int arcind=SOARCind[PH][thearc]; 
									sum1+=sz[sn][arcind];
								}
							}
						}
						int lastndind=SONDind[PH][thenode];
						model.add(su[sn][lastndind]-sv[sn][lastndind]+sum1+MINIdmd[pp][nn]<=TQ[thenode]);
						sum1.end();
					}
				}
			}
		}
		{   
			for (int p=PH+1;p!=P;++p){
				for (int pn=0;pn!= STGPATHN[PH+1][p];++pn){ 
					for (int sn=0;sn<STGSFSN[PH+1][p][pn]-1;++sn){
						int s1=STGSFSset[PH+1][p][pn][sn];
						for (int mm=sn+1;mm<STGSFSN[PH+1][p][pn];++mm){
							int s2=STGSFSset[PH+1][p][pn][mm];
							for (int an=0;an!=PAN[p];++an){  
								int thearc=PAset[p][an];
								int arcind=SOARCind[PH][thearc]; 
								model.add(sz[s1][arcind]==sz[s2][arcind]);
							}
						}
					}
				}
			}
		}
		{    
			IloExpr sum1(env);
			for (int sn=0;sn!=STTSN[PH+1];++sn){
				for (int nn=0;nn!=SONDN[PH];++nn){
					int thend=SONDset[PH][nn]; 
					sum1+=STtrpr[PH+1][sn]*COST1[thend]*su[sn][nn]; 
					sum1+=STtrpr[PH+1][sn]*COST2[thend]*sv[sn][nn]; 
				}
			}
			for (int sn=0;sn!=STTSN[PH+1];++sn){
				for (int nn=0;nn!=SOARCN[PH];++nn){
					int thearc=SOARCset[PH][nn]; 
					sum1+=STtrpr[PH+1][sn]*COST3[thearc]*sz[sn][nn]; 
				}
			}
			model.add(F>=sum1);
			sum1.end();
		}
		
		{
			IloCplex cplex(env);
			
			cplex.extract(model);
			cplex.setParam(IloCplex::Threads, 1);  
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			cplex.solve();
			
			double theobj=cplex.getObjValue();
			double fval=cplex.getValue(F);
			Ffval[PH][SCE]=fval;
			transcost=theobj-fval; 
			result=theobj;
			if (PH<P-1){
				for (int i=0;i!=I;++i){
					int theind=PNDindex[i][PTN[PH]-1]; 
					double theu=cplex.getValue(u[theind]);
					Inven[SCE][PH][i]=theu;
					double thev=cplex.getValue(v[theind]);
					UVval[SCE][PH][i]=thev;
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
					Arcflow[SCE][PH][nn]=thef; 
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

double SepCGP (int SCE, int PH, int TSC){
	double result=0;
	double thepro=PRO[PH][TSC]; // this is the probability
	IloEnv env; // this is the first stage env
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVarArray d22(env, C4N[PH], -IloInfinity, 0, ILOFLOAT); 
	IloNumVarArray d23(env, C5N[PH], -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray d24(env, C6N, -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray d25(env, C7N[PH], -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray d26(env, C8N, -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray d27(env, C9N[PH], -IloInfinity, 0, ILOFLOAT); 
	IloNumVarArray d28(env, C28N, -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray d29(env, C29N, -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray d30(env, C30N, -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray d31(env, C31N[PH], -IloInfinity, IloInfinity, ILOFLOAT); 
	IloNumVarArray d10(env, C10N[PH], -IloInfinity, 0, ILOFLOAT); 
	IloNumVarArray dbc(env, TRCN[PH], 0, IloInfinity, ILOFLOAT); 
	try{
		
		IloModel model(env);
		{   
			model.add(IloMaximize(env, obj)); 
		} 
		{   
			IloExpr Z23(env); 
			IloExpr Z24(env); 
			IloExpr Z25(env); 
			IloExpr Z26(env); 
			IloExpr Z27(env); 
			IloExpr Z28(env); 
			IloExpr Z29(env); 
			IloExpr Z30(env);
			IloExpr Z31(env);
			IloExpr Z10(env);  
			IloExpr BDZ(env); 
			int counter=0;
			for (int nn=1;nn!=PTN[PH];++nn){ 
				for (int i=0;i!=I1;++i){  
					
					Z23+=sup[i][PH][TSC][nn]*d23[counter]; 
					counter+=1;
				}
			}
			for (int i=0;i!=I1;++i){
				Z24+=sup[i][PH][TSC][0]*d24[i];
			}
			int counter25=0;
			for (int nn=1;nn!=PTN[PH];++nn){ 
				for (int i=0;i!=I2;++i){   
					
					Z25+=sup[i+I1][PH][TSC][nn]*d25[counter25];
					counter25+=1;
				}
			}
			for (int i=0;i!=I2;++i){
				
				Z26+=sup[i+I1][PH][TSC][0]*d26[i];	
			}
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				Z27+=TQ[thend]*d27[nn];
			}
			if (PH<P-1){  
				for (int nn=0;nn!=SPANN[PH];++nn){
					int thenode=SPANset[PH][nn];      
					Z10+=(TQ[thenode]-MINIdmd[PH][nn])*d10[nn];
				}
			}
			for (int bb=0;bb!=BN;++bb){
				Z28+=Bcap[bb]*d28[bb];
			}
			for (int i=0;i!=I;++i){
				if (PH==0){
					Z29+=iniv[i]*d29[i];
				}
				if (PH>0){
					Z29+=(Inven[SCE][PH-1][i])*d29[i];
				}
			}
			if (PH>0){
				for (int i=0;i!=I;++i){
					Z30+=(UVval[SCE][PH-1][i])*d30[i];
				}
			}
			if (PH>0){
				for (int nn=0;nn!=PSAN[PH-1];++nn){  
					
					Z31+=(Arcflow[SCE][PH-1][nn])*d31[nn];
				}
			}
			for (int nn=0;nn!=TRCN[PH];++nn){  
				BDZ+=TRcutval[PH][nn]*dbc[nn];
			}
			model.add(obj==Z23+Z24+Z25+Z26+Z27+Z10+Z28+Z29+Z30+Z31+BDZ);
			Z23.end();
			Z24.end();
			Z25.end();
			Z26.end();
			Z27.end();
			Z10.end();
			Z28.end();
			Z29.end();
			Z30.end();
			Z31.end();
			BDZ.end();
		}
		///// add constraints /////
		{   // Constraints for u variables:
			for (int nn=0;nn!=PNN[PH];++nn){
				IloExpr sum1(env);
				int thec23ind1=Nc5ind1[PH][nn];
				if (thec23ind1!=-1){
					sum1+=d23[thec23ind1];
				}
				int thec23ind2=Nc5ind2[PH][nn];
				if (thec23ind2!=-1){
					sum1+=-d23[thec23ind2];
				}
				int thec24ind=Nc6ind[PH][nn];
				if (thec24ind!=-1){
					sum1+=d24[thec24ind];
				}
				int thec25ind1=Nc7ind1[PH][nn];
				if (thec25ind1!=-1){
					sum1+=d25[thec25ind1];
				}
				int thec25ind2=Nc7ind2[PH][nn];
				if (thec25ind2!=-1){
					sum1+=-d25[thec25ind2];
				}
				int thec26ind=Nc8ind[PH][nn];
				if (thec26ind!=-1){
					sum1+=d26[thec26ind];
				}
				int thec27ind=Nc9ind[PH][nn];
				if (thec27ind!=-1){
					sum1+=d27[thec27ind];
				}
				if (PH<P-1){
					for (int ndn=0;ndn!=C10N[PH];++ndn){
						int thec10ind=Ndc10ind[PH][ndn][nn];
						if (thec10ind==ndn){ //if the nd is associated with the nn-th c10 constraint
							sum1+=d10[thec10ind];
						}
					}
				}
				int thenode=PNset[PH][nn];
				int thei=Nodeset[0][thenode];
				int thet=Nodeset[1][thenode];
				if (thet==PTset[PH][PTN[PH]-1]){  // if it is the last nd of a site in the period
					for (int bb=0;bb!=TRCN[PH];++bb){ // the benders cuts
						double thebce=TRucutset[PH][thei][bb]; // locate the coefficients with this variable in the primal p
						sum1+=-thebce*dbc[bb]; // associated with the nn-th node in the nd set    
					}
				}
				model.add(sum1<=COST1[thenode]);
				sum1.end();
			}
		}
		{   
			for (int nn=0;nn!=PNN[PH];++nn){
				IloExpr sum1(env);
				int thec23ind1=Nc5ind1[PH][nn];
				if (thec23ind1!=-1){
					sum1+=d23[thec23ind1];
				}
				int thec23ind2=Nc5ind2[PH][nn];
				if (thec23ind2!=-1){
					sum1+=-d23[thec23ind2];
				}
				int thec24ind=Nc6ind[PH][nn];
				if (thec24ind!=-1){
					sum1+=d24[thec24ind];
				}
				int thec25ind1=Nc7ind1[PH][nn];
				if (thec25ind1!=-1){
					sum1+=-d25[thec25ind1];
				}
				int thec25ind2=Nc7ind2[PH][nn];
				if (thec25ind2!=-1){
					sum1+=d25[thec25ind2];
				}
				int thec26ind=Nc8ind[PH][nn];
				if (thec26ind!=-1){
					sum1+=-d26[thec26ind];
				}
				if (PH<P-1){
					for (int ndn=0;ndn!=C10N[PH];++ndn){
						int thec10ind=Ndc10ind[PH][ndn][nn];
						if (thec10ind==ndn){ //if the nd is associated with the nn-th c10 constraint
							sum1+=-d10[thec10ind];
						}
					}
				}
				int thenode=PNset[PH][nn];
				int thei=Nodeset[0][thenode];
				int thet=Nodeset[1][thenode];
				if (thet==PTset[PH][PTN[PH]-1]){  // if it is the last nd of a site in the period
					for (int bb=0;bb!=TRCN[PH];++bb){ // the benders cuts
						double thebce=TRvcutset[PH][thei][bb]; // locate the coefficients with this variable in the primal p
						sum1+=-thebce*dbc[bb];     // associated with the nn-th node in the nd set
					}
				}
				model.add(sum1<=COST2[thenode]);
				sum1.end();
			}
		}
		{   // Constraints for z variables:
			for (int nn=0;nn!=PAN[PH];++nn){
				IloExpr sum1(env);
				int thec22ind=Ac4ind[PH][nn];
				if (thec22ind!=-1){ // if we have a bid for the arc
					sum1+=d22[thec22ind];
				}
				int thec23ind=Ac5ind[PH][nn];
				if (thec23ind!=-1){   // if the arc is associated with a c23
					sum1+=d23[thec23ind];
				}
				int thec24ind=Ac6ind[PH][nn];
				if (thec24ind!=-1){  // if the arc is associated with a c6
					sum1+=d24[thec24ind];
				}
				int thec25ind=Ac7ind[PH][nn];
				if (thec25ind!=-1){   // if the arc is associated with a c7
					sum1+=-d25[thec25ind];
				}
				int thec26ind=Ac8ind[PH][nn];
				if (thec26ind!=-1){  // if the arc is associated with a c8
					sum1+=-d26[thec26ind];
				}
				if (PH<P-1){
					for (int ndn=0;ndn!=C10N[PH];++ndn){
						int thec10ind=Arcc10ind[PH][ndn][nn];
						if (thec10ind==ndn){  //if the arc is associated with the nn-th c10 constraint
							sum1+=d10[thec10ind];
						}
					}
				}
				int thearc=PAset[PH][nn]; // locate the arc
				int theind=PSAindex[PH][thearc];
				if (theind!=-1){
					for (int bb=0;bb!=TRCN[PH];++bb){ // the benders cuts
						double thebce=TRzcutset[PH][theind][bb]; // locate the coefficients with this variable in the primal p
						sum1+=-thebce*dbc[bb]; // associated with the nnth arc in the arc set
					}
				}
				model.add(sum1<=COST3[thearc]);
				sum1.end();
			}
		}
		{   // constraints for ty vairiables 
			int counter=0; 
			for (int b=0;b!=BN;++b){ // enumerate ty variables 
				IloExpr sum1(env);
				for (int n=0;n!=PBAN[PH][b];++n){
					sum1+=-d22[counter];
					counter+=1;
				}
				for (int cc=0;cc!=TRCN[PH];++cc){ // the benders cuts
					double thebce=TRtycutset[PH][b][cc]; // locate the coefficients with this variable in the primal p
					sum1+=-thebce*dbc[cc];      // associated with the nn-th node in the nd set
				}
				sum1+=d28[b];
				model.add(sum1==0);
				sum1.end();
			}
		}
		{  // constriants for tu variables 
			for (int i=0;i!=I;++i){
				IloExpr sum1(env);
				int thec24ind=TNc6ind[PH][i];
				if (thec24ind!=-1){
					sum1+=-d24[thec24ind];
				}
				int thec26ind=TNc8ind[PH][i];
				if (thec26ind!=-1){
					sum1+=-d26[thec26ind];
				}
				int thec29ind=TNc29ind[PH][i];
				if (thec29ind!=-1){
					sum1+=d29[thec29ind];
				}
				model.add(sum1==0);
				sum1.end();
			}
		}
		{  // constriants for tv variables 
			for (int i=0;i!=I;++i){
				IloExpr sum1(env);
				int thec24ind=TNc6ind[PH][i];
				if (thec24ind!=-1){
					sum1+=-d24[thec24ind];
				}
				int thec26ind=TNc8ind[PH][i];
				if (thec26ind!=-1){
					sum1+=d26[thec26ind];
				}
				int thec30ind=TNc30ind[PH][i];
				if (thec30ind!=-1){
					sum1+=d30[thec30ind];
				}
				model.add(sum1==0);
				sum1.end();
			}
		}
		{ //constraints for tz variables
			if (PH>0){
				for (int zn=0;zn!=PSAN[PH-1];++zn){
					IloExpr sum1(env);
					int thec25ind=TAc7ind[PH][zn]; 
					if (thec25ind!=-1){
						sum1+=-d25[thec25ind];
					}
					int thec26ind=TAc8ind[PH][zn]; 
					if (thec26ind!=-1){
						sum1+=-d26[thec26ind];
					}
					int thec31ind=TAc31ind[PH][zn]; 
					if (thec31ind!=-1){
						sum1+=d31[thec31ind];
					}
					if (PH<P-1){
						for (int ndn=0;ndn!=C10N[PH];++ndn){
							int thec10ind=TArcc10ind[PH][ndn][zn];
							if (thec10ind==ndn){ 
								sum1+=d10[thec10ind];
							}
						}
					}
					int thearc=PSAset[PH-1][zn]; 
					int theind=PSAindex[PH][thearc]; 
					if (theind!=-1){
						for (int bb=0;bb!=TRCN[PH];++bb){ 
							double thebce=TRzcutset[PH][theind][bb]; 
							sum1+=-thebce*dbc[bb]; 
						}
					}
					model.add(sum1==0);
					sum1.end();
				}
			}
		}
		{ 
			IloExpr sum1(env);
			for (int bb=0;bb!=TRCN[PH];++bb){
				sum1+=dbc[bb];
			}
			model.add(sum1<=1);
			sum1.end();
		}
		{
			IloCplex cplex(env);
			cplex.extract(model);
			cplex.setParam(IloCplex::Threads,1);  
			cplex.setOut(env.getNullStream());
			cplex.setWarning(env.getNullStream());
			cplex.solve();
			result=cplex.getObjValue();
			int theph=PH-1; 
			if (theph==-1){ // this is the first period (we thus reach the biding phase)
				double thefix=result*thepro; 
				for (int bb=0;bb!=BN;++bb){
					double res=cplex.getValue(d28[bb]); 
					res=res*thepro; // this is the probability
					PBidcutset[bb]+=res; 
					thefix-=res*Bcap[bb];
				}
				PBidcutval+=thefix;
			}
			if (theph>=0){ 
				double thefix=result*thepro; 
				for (int bb=0;bb!=BN;++bb){
					double res=cplex.getValue(d28[bb]);
					res=res*thepro;
					PTRtycutset[bb]+=res; 
					thefix-=res*Bcap[bb];
				}
				for (int i=0;i!=I;++i){
					double res=cplex.getValue(d29[i]);
					res=res*thepro;
					PTRucutset[i]+=res;
					thefix-=res*(Inven[SCE][theph][i]);
				}
				for (int i=0;i!=I;++i){
					double res=cplex.getValue(d30[i]);
					res=res*thepro;
					PTRvcutset[i]+=res;
					thefix-=res*(UVval[SCE][theph][i]);
				}
				for (int nn=0;nn!=PSAN[theph];++nn){  
					double res=cplex.getValue(d31[nn]);
					res=res*thepro;
					PTRzcutset[theph][nn]+=res; 
					thefix-=res*(Arcflow[SCE][theph][nn]);
				}
				PTRcutval+=thefix;
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

double CGP (int SCE, int PH){ 
	double result=0;
	int theph=PH; 
	for (int bn=0;bn!=BN;++bn){
		PTRtycutset[bn]=0;}
	for (int in=0;in!=I;++in){
		PTRucutset[in]=0;}
	for (int in=0;in!=I;++in){
		PTRvcutset[in]=0;}
	for (int nn=0;nn!=PSAN[theph];++nn){   
		PTRzcutset[theph][nn]=0; 
	}
	PTRcutval=0;
	for (int sn=0;sn!=SN;++sn){
		double theres=SepCGP (SCE, PH+1, sn);
		result+=theres*PRO[PH+1][sn];
	}
	TRCN[theph]+=1;
	for (int bb=0;bb!=BN;++bb){
		TRtycutset[theph][bb].resize(TRCN[theph],0); 
	}
	for (int i=0;i!=I;++i){
		TRucutset[theph][i].resize(TRCN[theph],0);
	}
	for (int i=0;i!=I;++i){
		TRvcutset[theph][i].resize(TRCN[theph],0);
	}
	for (int nn=0;nn!=PSAN[theph];++nn){  
		TRzcutset[theph][nn].resize(TRCN[theph],0); 
	}
	TRcutval[theph].resize(TRCN[theph],0);
	int thecut=TRCN[theph]-1;
	for (int bb=0;bb!=BN;++bb){
		TRtycutset[theph][bb][thecut]+=PTRtycutset[bb]; 
	}
	for (int i=0;i!=I;++i){
		TRucutset[theph][i][thecut]+=PTRucutset[i];
	}
	for (int i=0;i!=I;++i){
		TRvcutset[theph][i][thecut]+=PTRvcutset[i];
	}
	for (int nn=0;nn!=PSAN[theph];++nn){  
		TRzcutset[theph][nn][thecut]+=PTRzcutset[theph][nn]; 
	}
	TRcutval[theph][thecut]+=PTRcutval;
	return result;
}

double CGPS0(){
	double btemp=0;
	double Presult=0;
	for (int bn=0;bn!=BN;++bn){
		PBidcutset[bn]=0;}
	PBidcutval=0;
	Presult=0;
	for (int sn=0;sn!=SN;++sn){
		double theres=SepCGP (0,0, sn);
		Presult+=theres*PRO[0][sn];
	}
	BidCN+=1; 
	for (int bb=0;bb!=BN;++bb){
		Bidcutset[bb].resize(BidCN,0); 
	}
	Bidcutval.resize(BidCN,0); 
	int thecut=BidCN-1;
	for (int bb=0;bb!=BN;++bb){
		Bidcutset[bb][thecut]+=PBidcutset[bb];
	}
	Bidcutval[thecut]+=PBidcutval;
	btemp+=Presult;
	return btemp;
}

void Backward(int isint){ 
	Bgap=pdlgap; 
	double btemp=-1;
	double ftemp=-1;
	for (int ss=0;ss!=TSN;++ss){  
		for (int pp=P-2;pp>=0;--pp){  
			double crtgap=100;
			ftemp=Ffval[pp][ss];
			btemp=CGP(ss,pp);    
			crtgap=(btemp-ftemp)/btemp*double(100);
			while (crtgap>Bgap){
				TRP (ss,pp);         // solve the updated problem
				ftemp=Ffval[pp][ss];
				btemp=CGP(ss,pp);    //backward problem (solve the next stage problems)
				crtgap=(btemp-ftemp)/btemp*double(100);
			}
		}
	}
	if (isint==0){Bgap=pdlgap;}
	if (isint==1){Bgap=pdlgap02;}
	double crtgap=100;
	ftemp=Bidf;
	btemp=CGPS0();
	crtgap=(btemp-ftemp)/btemp*double(100);
	cout<<"Starting Stage-0 Bgap=["<<Bgap<<"] crtgap=["<<crtgap<<"]"<<endl;
	if (isint==0){
		while (crtgap>Bgap){
			BIDP (1,1);
			ftemp=Bidf;
			btemp=CGPS0();
			crtgap=(btemp-ftemp)/btemp*double(100);
		}
	}
	if (isint==1){
		while (crtgap>Bgap){
			BIDP (1,2);
			ftemp=Bidf;
			btemp=CGPS0();
			crtgap=(btemp-ftemp)/btemp*double(100);
		}
	}	
}

void Forward(){ 
	for (int ss=0;ss!=TSN;++ss){
		for (int pp=0;pp!=P-1;++pp){
			TRP(ss,pp);
		}
	}
}

double SDDP(double err,int stn,int isint){
	int stable=0;    
	double tempgap; 
	double OLB=0;    
	double theresult; 
	double flb=0;     
	struct timeval t1;   
	struct timeval t2; 
	double uptime;
	gettimeofday(&t1,NULL);
	if (isint==0){
		uptime=uptime_p1;
	}
	if (isint==1){
		uptime=uptime_p2;
	}
	int itrcounter=0;
	while (stable<stn){
		itrcounter+=1;
		Sampling(TSN);
		if (isint==0){
			theresult=BIDP(1,0);
		    Forward(); 
			Backward(isint); 
		}
		if (isint==1){
			theresult= FLOATBIDP(1,timlim01);
			Forward(); 
			Backward(isint); 
		} 
		int addstep=1;
		theresult=thebidlb; 
		if (LB<theresult){
			LB=theresult;
			flb=LB;
		}
		if (OLB==0){
			tempgap=1;}
		else{
			tempgap=double(theresult-OLB)/double(OLB);}
		tempgap=tempgap*double(100);
		cout<<"SDDPthegap="<<tempgap<<endl;
		OLB=theresult; 
		if (tempgap<err){  
			stable+=addstep;}
		else{
			stable=0;}
		gettimeofday(&t2,NULL);
		if ((((t2.tv_sec - t1.tv_sec)*1000000+(t2.tv_usec - t1.tv_usec))/double(1000000))>=uptime){
			cout<<"time up"<<endl;
			break;
		}		
	}
	if (isint==0){
        itrn_p1=itrcounter;
	}
	if (isint==1){
        itrn_p2=itrcounter;
	}
	cout<<"thelb="<<flb<<endl;
	cout<<endl;
	cout<<endl;
	cout<<endl;
	return flb; 
}

double TRP02 (int PH){   
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
	IloNumVar F(env, 0, IloInfinity, ILOFLOAT); 
	try{
		IloModel model(env);
		{   
			model.add(IloMinimize(env, obj)); 
		} 
		{   
			IloExpr Z2(env); 
			IloExpr Z3(env); 
			IloExpr Z4(env); 
			IloExpr Z5(env); 
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
			model.add(obj==Z2+Z3+Z4+F);
			Z2.end();
			Z3.end();
			Z4.end();
			Z5.end();
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
					
					model.add(u[ndind1]+v[ndind1]==u[ndind2]+v[ndind2]+PD[nd1]-sum1); 
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
				model.add(u[ndind1]+v[ndind1]==tu[i]+tv[i]+PD[nd1]-sum1);  
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
					model.add(u[ndind1]-v[ndind1]==u[ndind2]-v[ndind2]+PD[nd1]+tsum1+sum1); 
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
				model.add(u[ndind1]-v[ndind1]==tu[i+I1]-tv[i+I1]+PD[nd1]+tsum1+sum1); 
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
					model.add(tu[i]==PInven[PH-1][i]); 
				}
			}
		}
		{   
			for (int i=0;i!=I;++i){
				if (PH==0){
					model.add(tv[i]==0); 
				}
				if (PH>0){
					model.add(tv[i]==PUVval[PH-1][i]); 
				}
			}
		}
		{   
			if (PH>0){
				for (int nn=0;nn!=PSAN[PH-1];++nn){  
					model.add(tz[nn]==PArcflow[PH-1][nn]);
				}
			}
		}
		{   
			if (PH<P-1){ 
				for (int n=0;n!=TRCN[PH];++n){ 
					IloExpr sum1(env);
					for (int i=0;i!=I;++i){
						double theval=TRucutset[PH][i][n]; 
						int un=PNDindex[i][PTN[PH]-1];     
						sum1+=theval*u[un];
					}
					for (int i=0;i!=I;++i){
						double theval=TRvcutset[PH][i][n]; 
						int un=PNDindex[i][PTN[PH]-1];     
						sum1+=theval*v[un];
					}
					for (int zn=0;zn!=PSAN[PH];++zn){
						double theval=TRzcutset[PH][zn][n]; 
						int thearc=PSAset[PH][zn];   
						int thep=ARCPN[thearc]; 
						if (thep==PH){  
							int theind=ARCPS[thearc];  
							sum1+=theval*z[theind];
						}
						else {    
							int theind=PSAindex[PH-1][thearc];
							sum1+=theval*tz[theind]; 
						}
					}	
					for (int yn=0;yn!=BN;++yn){
						double theval=TRtycutset[PH][yn][n]; 
						sum1+=theval*ty[yn];
					}
					model.add(F>=sum1+TRcutval[PH][n]);
					sum1.end();
				}
			}
		}
		
		if (PH<P-1){ 
			IloArray<IloNumVarArray> sz(env, STTSN[PH+1]);
			IloArray<IloNumVarArray> su(env, STTSN[PH+1]);
			IloArray<IloNumVarArray> sv(env, STTSN[PH+1]);
			for (int sn=0;sn!=STTSN[PH+1];++sn){ 
				sz[sn]=IloNumVarArray (env, SOARCN[PH], 0, IloInfinity, ILOFLOAT); 
				su[sn]=IloNumVarArray (env, SONDN[PH], 0, IloInfinity, ILOFLOAT); 
				sv[sn]=IloNumVarArray (env, SONDN[PH], 0, IloInfinity, ILOFLOAT); 
			}
			{   
				for (int sn=0;sn!=STTSN[PH+1];++sn){ 
					for (int b=0;b!=BN;++b){
						for (int pp=PH+1;pp!=P;++pp){  
							for (int n=0;n!=PBAN[pp][b];++n){
								int thearc=PBAset[pp][b][n];
								int theind=SOARCind[PH][thearc]; 
								model.add(sz[sn][theind]<=ty[b]);
							}
						}
					}
				}
			}
			{   
				for (int sn=0;sn!=STTSN[PH+1];++sn){
					for (int t1=PTset[PH+1][1];t1!=T;++t1){  
						int t2=t1-1; 
						for (int i=0;i!=I1;++i){  
							int nd1=NDindex[i][t1]; 
							int nd2=NDindex[i][t2];
							int ndind1=SONDind[PH][nd1];  
							int ndind2=SONDind[PH][nd2];  
							IloExpr sum1(env);
							for (int n=0;n!=OUTAN[nd1];++n){
								int thearc=OUTAset[nd1][n];
								int theind=SOARCind[PH][thearc]; 
								sum1+=sz[sn][theind];
							}
							model.add(su[sn][ndind1]+sv[sn][ndind1]==su[sn][ndind2]+sv[sn][ndind2]+STD[PH+1][sn][nd1]-sum1); 
							sum1.end();
						}
					}
				}			
			}
			{   
				for (int sn=0;sn!=STTSN[PH+1];++sn){
					int t1=PTset[PH+1][0];  
					for (int i=0;i!=I1;++i){
						int oldndind=PNDindex[i][PTN[PH]-1]; 
						int nd1=NDindex[i][t1]; 
						int ndind=SONDind[PH][nd1]; 
						IloExpr sum1(env);
						for (int n=0;n!=OUTAN[nd1];++n){
							int thearc=OUTAset[nd1][n];
							int theind=SOARCind[PH][thearc]; 
							sum1+=sz[sn][theind];
						}
						
						model.add(su[sn][ndind]+sv[sn][ndind]==u[oldndind]+v[oldndind]+STD[PH+1][sn][nd1]-sum1);  
						sum1.end();
					}
				}
			}
			{   
				for (int sn=0;sn!=STTSN[PH+1];++sn){
					for (int t1=PTset[PH+1][1];t1!=T;++t1){ 
						int t2=t1-1;
						for (int i=0;i!=I2;++i){   
							int nd1=NDindex[i+I1][t1]; 
							int nd2=NDindex[i+I1][t2];
							IloExpr sum1(env);  
							for (int n=0;n!=INAN[nd1];++n){
								int thearc=INAset[nd1][n];
								int thephase=ARCPN[thearc]; 
								if (thephase<PH){ 
									int theind=PSAindex[PH-1][thearc]; 
									sum1+=tz[theind];
								}
								if (thephase==PH){ 
									int theind=ARCPS[thearc]; 
									sum1+=z[theind];
								}
								if (thephase>PH){ 
									int arcind=SOARCind[PH][thearc]; 
									sum1+=sz[sn][arcind];
								}
							}
							int ndind1=SONDind[PH][nd1];  
							int ndind2=SONDind[PH][nd2];  
							model.add(su[sn][ndind1]-sv[sn][ndind1]==su[sn][ndind2]-sv[sn][ndind2]+STD[PH+1][sn][nd1]+sum1); 
							sum1.end();
						}
					}
				}
			}
			{   
				for (int sn=0;sn!=STTSN[PH+1];++sn){
					int t1=PTset[PH+1][0];  
					for (int i=0;i!=I2;++i){   
						int oldndind=PNDindex[i+I1][PTN[PH]-1]; 
						int nd1=NDindex[i+I1][t1]; 
						IloExpr sum1(env);  
						for (int n=0;n!=INAN[nd1];++n){
							int thearc=INAset[nd1][n];
							int thephase=ARCPN[thearc]; 
							if (thephase<PH){ 
								int theind=PSAindex[PH-1][thearc]; 
								sum1+=tz[theind];
							}
							if (thephase==PH){ 
								int theind=ARCPS[thearc]; 
								sum1+=z[theind];
							}
							if (thephase>PH){ 
								int arcind=SOARCind[PH][thearc]; 
								sum1+=sz[sn][arcind];
							}
						}
						int ndind1=SONDind[PH][nd1];  
						model.add(su[sn][ndind1]-sv[sn][ndind1]==u[oldndind]-v[oldndind]+STD[PH+1][sn][nd1]+sum1); 
						sum1.end();
					}
				}
			}
			{   
				for (int sn=0;sn!=STTSN[PH+1];++sn){
					for (int nn=0;nn!=SONDN[PH];++nn){
						int thend=SONDset[PH][nn]; 
						model.add(su[sn][nn]<=TQ[thend]);
					}
				}
			}
			{ 
				for (int sn=0;sn!=STTSN[PH+1];++sn){
					for (int pp=PH+1;pp!=P-1;++pp){
						for (int nn=0;nn!=SPANN[pp];++nn){
							
							int thenode=SPANset[pp][nn];      
							int thesite=Nodeset[0][thenode];  
							int lastnd=NDindex[thesite][PTN[pp]-1]; 
							int thetime=Nodeset[1][thenode];  
							IloExpr sum1(env);
							for (int tt=PTset[pp+1][0]; tt<=thetime;++tt){
								int crtnd=NDindex[thesite][tt]; 
								for (int n=0;n!=INAN[crtnd];++n){
									int thearc=INAset[crtnd][n];
									int thephase=ARCPN[thearc]; 
									if (thephase<PH){ 
										int theind=PSAindex[PH-1][thearc]; 
										sum1+=tz[theind];
									}
									if (thephase==PH){ 
										int theind=ARCPS[thearc]; 
										sum1+=z[theind];
									}
									if (thephase>PH){ 
										int arcind=SOARCind[PH][thearc]; 
										sum1+=sz[sn][arcind];
									}
								}
							}
							int lastndind=SONDind[PH][thenode];
							model.add(su[sn][lastndind]-sv[sn][lastndind]+sum1+MINIdmd[pp][nn]<=TQ[thenode]);
							sum1.end();
						}
					}
				}
			}
			{   
				for (int p=PH+1;p!=P;++p){
					for (int pn=0;pn!= STGPATHN[PH+1][p];++pn){ 
						for (int sn=0;sn<STGSFSN[PH+1][p][pn]-1;++sn){
							int s1=STGSFSset[PH+1][p][pn][sn];
							for (int mm=sn+1;mm<STGSFSN[PH+1][p][pn];++mm){
								int s2=STGSFSset[PH+1][p][pn][mm];
								for (int an=0;an!=PAN[p];++an){  
									int thearc=PAset[p][an];
									int arcind=SOARCind[PH][thearc]; 
									model.add(sz[s1][arcind]==sz[s2][arcind]);
								}
							}
						}
					}
				}
			}
			{    
				IloExpr sum1(env);
				for (int sn=0;sn!=STTSN[PH+1];++sn){
					for (int nn=0;nn!=SONDN[PH];++nn){
						int thend=SONDset[PH][nn]; 
						sum1+=STtrpr[PH+1][sn]*COST1[thend]*su[sn][nn]; 
						sum1+=STtrpr[PH+1][sn]*COST2[thend]*sv[sn][nn]; 
					}
				}
				for (int sn=0;sn!=STTSN[PH+1];++sn){
					for (int nn=0;nn!=SOARCN[PH];++nn){
						int thearc=SOARCset[PH][nn]; 
						sum1+=STtrpr[PH+1][sn]*COST3[thearc]*sz[sn][nn]; 
					}
				}
				model.add(F>=sum1);
				sum1.end();
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
			double fval=cplex.getValue(F);
			transcost=theobj-fval; 
			result=transcost; 
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				double theval=cplex.getValue(u[nn]);
				TTivvol+=theval;
				if(Nodeset[0][thend]<I1){
					TTivvol_s+=theval;
				}
				else{
				    TTivvol_d+=theval;
				}
				TTivcost+=theval*COST1[thend]; 
			}
			for (int nn=0;nn!=PNN[PH];++nn){
				int thend=PNset[PH][nn]; 
				double theval=cplex.getValue(v[nn]);
				TTbgvol+=theval;
				if(Nodeset[0][thend]<I1){
					TTbgvol_s+=theval;
				}
				else{
					TTbgvol_d+=theval;
				}
				TTbgcost+=theval*COST2[thend]; 
			}
			for (int nn=0;nn!=PAN[PH];++nn){
				int thearc=PAset[PH][nn];
				double theval=cplex.getValue(z[nn]);
				if (COST3[thearc]<0.001){
					TTbidvol+=theval;
				}
				if (COST3[thearc]>0.001){
					TTspotvol+=theval;
					TTspotcost+=theval*COST3[thearc];
				}
			}
			if (PH<P-1){
				for (int i=0;i!=I;++i){
					int theind=PNDindex[i][PTN[PH]-1]; 
					double theu=cplex.getValue(u[theind]);
					PInven[PH][i]=theu;
					double thev=cplex.getValue(v[theind]);
					PUVval[PH][i]=thev;
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
					PArcflow[PH][nn]=thef; 
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

double SCECAl(int SIND){
	for (int p=0;p!=P;++p){
        int thess=random(SN); // randomly generate the sample scenario
		for (int n=0;n!=PTN[p];++n){  
			int thetime=PTset[p][n];
			for (int i=0;i!=I;++i){
				int thend=NDindex[i][thetime];
				PD[thend]=sup[i][p][thess][n]; 
			}
		}
	}
	double theupper=bidcost;
	double scetran=0;
	for (int p=0;p!=P-1;++p){ 
		for (int i=0;i!=I;++i){
			PInven[p][i]=0;
			PUVval[p][i]=0;
		}
		PArcflow[p].resize(PSAN[p],0);
	}
	for (int pp=0;pp!=P;++pp){
		double transcost=TRP02(pp);
		scetran+=transcost; 
	}
	theupper+=scetran;
	return theupper;
}

void Upperbound(int M){ 
	vector<double> Fres(M);
	UB=0;
	int STSN=0;
	for (int mm=0;mm!=M;++mm){
		double theres=SCECAl(mm);
		UB+=theres;
		Fres[STSN]=theres;
		STSN+=1;
	}
	double aveub=double(UB)/double(M);
	double totaldev=0;
	for (int nn=0;nn!=M;++nn){
		totaldev+=(double(Fres[nn]-aveub))*(double(Fres[nn]-aveub));
	}
	double avedev=(totaldev)/(M-1);
	UB=aveub+1.96*pow(avedev,0.5)/pow(double(M),0.5);
	stgap=double(UB-LB)/double(UB);
	stgap*=100;
	for (int b=0;b!=BN;++b){
		TTbidcap+=Bcap[b]*SHN[b];
	}
	TTspotcost=TTspotcost/double(M);
	TTivcost=TTivcost/double(M);
	TTbgcost=TTbgcost/double(M);
	TTbidvol=TTbidvol/double(M);
	TTspotvol=TTspotvol/double(M);
	TTivvol=TTivvol/double(M);
	TTbgvol=TTbgvol/double(M);
	TTbirate=TTbirate/double(M);
	TTivvol_s=TTivvol_s/double(M);
	TTivvol_d=TTivvol_d/double(M);
	TTbgvol_s=TTbgvol_s/double(M);
	TTbgvol_d=TTbgvol_d/double(M);
}

void Enumeration() {
	cout << "bidcost=" << bidcost << endl;
	int search[P];
	for (int p = 0; p != P - 1; ++p) {
		search[p] = 0;
	}
	UB = 0;
	int STSN = 0;
	int doable = 1;
	while (doable == 1) {
		int thes = 0;
		while (thes < SN) {
			search[P - 1] = thes; // load the last stage scenario
			for (int p = 0; p != P - 1; ++p) { // this circulation does not include the last stage 
				for (int n = 0; n != PTN[p]; ++n) {  // record the supply/demand
					int thetime = PTset[p][n];
					for (int i = 0; i != I; ++i) {
						int thend = NDindex[i][thetime];// thend
						PD[thend] = sup[i][p][search[p]][n];  // we solve scenario by scenario, each time only one scenario is loaded into the tree (as the first scenario)
					}
				}
			}
			for (int n = 0; n != PTN[P - 1]; ++n) {   // record the supply/demand of last stage
				int thetime = PTset[P - 1][n];
				for (int i = 0; i != I; ++i) {
					int thend = NDindex[i][thetime];// thend
					PD[thend] = sup[i][P - 1][thes][n];
				}
			}
			for (int p = 0; p != P - 1; ++p) {
				for (int i = 0; i != I; ++i) {
					PInven[p][i] = 0;
					PUVval[p][i] = 0;
				}
				PArcflow[p].resize(PSAN[p], 0);
			}
			STSN += 1;
			double theupper = bidcost;
			double scetran = 0;
			for (int pp = 0; pp != P; ++pp) {
				double transcost = TRP02(pp);
				scetran += transcost;
			}
			theupper += scetran;
			UB += theupper;
			thes += 1;
			if (thes >= SN) {
				break;
			}
		}
		doable = 0;
		for (int p = P - 2; p >= 0; --p) {
			if (search[p] < SN - 1) {
				search[p] += 1;
				doable = 1;
				for (int pp = p + 1; pp < P - 1; ++pp) {
					search[pp] = 0;
				}
				break;
			}
		}
	}
	cout << "STSN=" << STSN << endl;
	double aveub = double(UB) / double(STSN);
	UB = aveub;
	cout << "aveub para=" << aveub << endl;
	stgap=double(aveub-LB)/double(aveub);
	stgap *= 100;
	for (int b = 0; b != BN; ++b) {
		TTbidcap += Bcap[b] * SHN[b];
	}
	TTspotcost = TTspotcost / double(STSN);
	TTivcost = TTivcost / double(STSN);
	TTbgcost = TTbgcost / double(STSN);
	TTbidvol = TTbidvol / double(STSN);
	TTspotvol = TTspotvol / double(STSN);
	TTivvol = TTivvol / double(STSN);
	TTbgvol = TTbgvol / double(STSN);
	TTbirate = TTbirate / double(STSN);
	TTivvol_s = TTivvol_s / double(STSN);
	TTivvol_d = TTivvol_d / double(STSN);
	TTbgvol_s = TTbgvol_s / double(STSN);
	TTbgvol_d = TTbgvol_d / double(STSN);
}

void ORISDDP(double err){
    struct timeval p1_start,p2_start,p3_start,p4_start;   
	struct timeval p1_end,p2_end,p3_end,p4_end;   
	gettimeofday(&p1_start,NULL);
	lb_p1=SDDP(err,10,0); 
	gettimeofday(&p1_end,NULL);
    t_p1=((p1_end.tv_sec - p1_start.tv_sec)*1000000+(p1_end.tv_usec - p1_start.tv_usec))/double(1000000);
    tbcn_p1=BidCN;
	for (int nn=0;nn!=P-1;++nn){
		tbcn_p1+=TRCN[nn];
	}
	gettimeofday(&p2_start,NULL);
	//lb_p2=SDDP(err,10,1); 
	tbcn_total=BidCN;
	for (int nn=0;nn!=P-1;++nn){
		tbcn_total+=TRCN[nn];
	}
	tbcn_p2=tbcn_total-tbcn_p1;
	gettimeofday(&p2_end,NULL);
    t_p2=((p2_end.tv_sec - p2_start.tv_sec)*1000000+(p2_end.tv_usec - p2_start.tv_usec))/double(1000000);
	gettimeofday(&p3_start,NULL);
	LB= FLOATBIDP(1,timlim02); 
	gettimeofday(&p3_end,NULL);
    t_p3=((p3_end.tv_sec - p3_start.tv_sec)*1000000+(p3_end.tv_usec - p3_start.tv_usec))/double(1000000);
    t_total=((p3_end.tv_sec - p1_start.tv_sec)*1000000+(p3_end.tv_usec - p1_start.tv_usec))/double(1000000);
	lb_p3=LB;
	lb_total=LB;
	cout<<"LB="<<LB<<endl;
	if (P<=7){
		ALLSN=pow(double(SN),double(P));}
	else {
		ALLSN=999999;}
	gettimeofday(&p4_start,NULL);
	if (ALLSN<=ESTV){ 
	    Enumeration();}
	if (ALLSN>ESTV){
		Upperbound(ESTV);
	}
	gettimeofday(&p4_end,NULL);
	t_p4=((p4_end.tv_sec - p4_start.tv_sec)*1000000+(p4_end.tv_usec - p4_start.tv_usec))/double(1000000);
	itrn_total=itrn_p1+itrn_p2;
	upperbound=UB;
	totalgap=stgap;
}

void documentation(){
	cout <<"succeeded" << endl;
	if (ALLSN<=ESTV){ 
		cout<<"Real UB"<<endl;}
	if (ALLSN>ESTV){
		cout<<"Statistical UB"<<endl;
	}
	cout<<"upperbound="<<upperbound<<endl;
	cout<<"totalgap="<<totalgap<<endl;
	cout<<"lb_total="<<lb_total<<endl;
	cout<<"lb_p1="<<lb_p1<<endl;
	cout<<"lb_p2="<<lb_p2<<endl;
	cout<<"lb_p3="<<lb_p3<<endl;
	cout<<"t_total="<<t_total<<endl;
	cout<<"t_p1="<<t_p1<<endl;
	cout<<"t_p2="<<t_p2<<endl;
	cout<<"t_p3="<<t_p3<<endl;
	cout<<"t_p4="<<t_p4<<endl;
	cout<<"itrn_total="<<itrn_total<<endl;
	cout<<"itrn_p1="<<itrn_p1<<endl;
	cout<<"itrn_p2="<<itrn_p2<<endl;
	cout<<"tbcn_total="<<tbcn_total<<endl;
	cout<<"tbcn_p1="<<tbcn_p1<<endl;
	cout<<"tbcn_p2="<<tbcn_p2<<endl;
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
	cout<<"Tivvol_s="<<TTivvol_s<<endl;
	cout<<"Tivvol_d="<<TTivvol_d<<endl;
	cout<<"Tbgvol_s="<<TTbgvol_s<<endl;
	cout<<"Tbgvol_d="<<TTbgvol_d<<endl;

	
	if (ALLSN<=ESTV){ 
		mycout01<<"Real UB"<<endl;}
	if (ALLSN>ESTV){
		mycout01<<"Statistical UB"<<endl;
	}
	mycout01<<"upperbound="<<endl<<upperbound<<endl;
	mycout01<<"totalgap="<<endl<<totalgap<<endl;
	mycout01<<"lb_total="<<endl<<lb_total<<endl;
	mycout01<<"t_total="<<endl<<t_total<<endl;
	mycout01<<"lb_p1="<<lb_p1<<endl;
	mycout01<<"lb_p2="<<lb_p2<<endl;
	mycout01<<"lb_p3="<<lb_p3<<endl;
	mycout01<<"t_p1="<<t_p1<<endl;
	mycout01<<"t_p2="<<t_p2<<endl;
	mycout01<<"t_p3="<<t_p3<<endl;
	mycout01<<"t_p4="<<t_p4<<endl;
	mycout01<<"itrn_total="<<itrn_total<<endl;
	mycout01<<"itrn_p1="<<itrn_p1<<endl;
	mycout01<<"itrn_p2="<<itrn_p2<<endl;
	mycout01<<"tbcn_total="<<tbcn_total<<endl;
	mycout01<<"tbcn_p1="<<tbcn_p1<<endl;
	mycout01<<"tbcn_p2="<<tbcn_p2<<endl;
	mycout01<<"bidcost="<<bidcost<<endl;
	for (int bb=0;bb!=BN;++bb){
		if (Bcap[bb]>1){
			mycout01<<"bid["<<bb<<"] purchased at"<<Bcap[bb]<<endl;
		}
	}
	mycout01<<"Tbidcap="<<TTbidcap<<endl;
	mycout01<<"Tbidcost="<<bidcost<<endl;
	mycout01<<"Tspotcost="<<TTspotcost<<endl;
	mycout01<<"Tivcost="<<TTivcost<<endl;
	mycout01<<"Tbgcost="<<TTbgcost<<endl;
	mycout01<<"Tbidvol="<<TTbidvol<<endl;
	mycout01<<"Tspotvol="<<TTspotvol<<endl;
	mycout01<<"Tivvol="<<TTivvol<<endl;
	mycout01<<"Tbgvol="<<TTbgvol<<endl;
	mycout01<<"Tivvol_s="<<TTivvol_s<<endl;
	mycout01<<"Tivvol_d="<<TTivvol_d<<endl;
	mycout01<<"Tbgvol_s="<<TTbgvol_s<<endl;
	mycout01<<"Tbgvol_d="<<TTbgvol_d<<endl;
}

int main (int argc, char **argv) {
	srand(2);
	srand((int)time(NULL));
	inputData();
	cout << "formulation" << endl;
	ORISDDP(sddperr);
    documentation();
	return 0;
}
