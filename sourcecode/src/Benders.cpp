#include <ilcplex/ilocplex.h>
#include <stdio.h>            
#include <stdlib.h>      
#include <sys/time.h>              /* time */         
#include <fstream>  
#include <string>
#include <math.h>
#include <sstream>
#include <algorithm>     
#include <vector>
#define random(x) (rand()%x)
namespace patch { template < typename T > std::string to_string(const T& n) { std::ostringstream stm; stm << n; return stm.str(); } }

using namespace std;
//input data starts here
int Dev = 10;// Deviation
extern const int CNN = 1;// Case Index
extern const int P = 3;// Number of phases
extern const int PT = 6;// Number of periods in a phase
extern const int T = 18;// Number of periods
extern const int I1 = 1;// Number of suppliers
extern const int I2 = 8;// Number of customers
extern const int I = 9;// Number of sites
extern const int BN = 72;//Number of bids
extern const int LBN = 9;//Number of bids on each arc
extern const int SPN = 339;//Number of shipments
extern const int SN = 10;//Number of scenarios at each stage
extern const int MBSN = 9;//Number of shipments in a bid
// tree parameters:
double PRO[P][SN] = { {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1},{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1},{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1} }; // the probobility of each scenario at each phase 
// timing parameters: 
int PTN[P] = { 6,6,6 }; // number of periods in each phase 
int PTset[P][PT] = { {0,1,2,3,4,5},{6,7,8,9,10,11},{12,13,14,15,16,17} }; // set of periods in each phase 
// distribution parameters: 
int sup[I][P][SN][PT] = { {{{11408,0,0,0,0,0},{11615,0,0,0,0,0},{11400,0,0,0,0,0},{11421,0,0,0,0,0},{11558,0,0,0,0,0},{11534,0,0,0,0,0},{11594,0,0,0,0,0},{11595,0,0,0,0,0},{11726,0,0,0,0,0},{11553,0,0,0,0,0}},{{11316,0,0,0,0,0},{11677,0,0,0,0,0},{11516,0,0,0,0,0},{11630,0,0,0,0,0},{11392,0,0,0,0,0},{11482,0,0,0,0,0},{11441,0,0,0,0,0},{11563,0,0,0,0,0},{11761,0,0,0,0,0},{11582,0,0,0,0,0}},{{11667,0,0,0,0,0},{11389,0,0,0,0,0},{11415,0,0,0,0,0},{11419,0,0,0,0,0},{11460,0,0,0,0,0},{11422,0,0,0,0,0},{11466,0,0,0,0,0},{11689,0,0,0,0,0},{11472,0,0,0,0,0},{11566,0,0,0,0,0}}},{{{-185,-221,-193,-187,-199,-209},{-195,-197,-201,-185,-207,-209},{-189,-183,-185,-209,-219,-201},{-199,-183,-185,-203,-209,-219},{-195,-195,-223,-193,-217,-201},{-221,-195,-207,-197,-213,-201},{-203,-217,-215,-183,-215,-187},{-217,-211,-197,-217,-213,-187},{-213,-223,-189,-207,-221,-217},{-215,-215,-185,-201,-183,-189}},{{-213,-191,-199,-183,-185,-215},{-195,-219,-213,-203,-209,-215},{-195,-215,-209,-219,-221,-221},{-205,-189,-223,-189,-189,-205},{-203,-199,-201,-199,-187,-215},{-203,-205,-201,-199,-215,-213},{-203,-191,-191,-205,-207,-183},{-193,-187,-223,-213,-197,-223},{-223,-199,-223,-197,-191,-223},{-191,-203,-195,-193,-219,-217}},{{-215,-223,-217,-221,-195,-205},{-203,-215,-199,-193,-187,-213},{-209,-207,-221,-207,-183,-189},{-215,-221,-189,-217,-191,-207},{-195,-189,-197,-223,-199,-205},{-191,-207,-203,-207,-211,-187},{-223,-211,-215,-193,-183,-195},{-211,-201,-213,-221,-207,-217},{-207,-201,-223,-197,-219,-191},{-189,-205,-209,-205,-223,-195}}},{{{-172,-156,-148,-169,-163,-150},{-161,-167,-171,-159,-167,-167},{-147,-147,-174,-169,-155,-156},{-159,-161,-155,-163,-159,-150},{-161,-169,-148,-177,-153,-175},{-145,-172,-177,-156,-156,-163},{-153,-177,-145,-166,-177,-169},{-150,-164,-159,-150,-156,-172},{-158,-172,-167,-159,-164,-167},{-147,-150,-159,-177,-171,-172}},{{-150,-164,-167,-155,-171,-161},{-156,-175,-163,-163,-177,-145},{-148,-155,-167,-166,-155,-155},{-175,-145,-167,-175,-145,-150},{-153,-171,-169,-159,-161,-172},{-172,-147,-156,-158,-156,-172},{-147,-150,-148,-147,-164,-163},{-174,-167,-147,-163,-148,-164},{-158,-166,-158,-167,-175,-145},{-172,-169,-150,-156,-172,-174}},{{-164,-158,-148,-156,-171,-158},{-150,-172,-156,-172,-145,-172},{-151,-164,-163,-175,-169,-161},{-161,-161,-158,-151,-163,-163},{-174,-159,-161,-148,-164,-156},{-167,-153,-159,-163,-145,-147},{-155,-169,-158,-153,-175,-163},{-158,-169,-151,-175,-171,-163},{-175,-156,-150,-174,-167,-167},{-166,-163,-155,-159,-153,-150}}},{{{-231,-229,-243,-208,-238,-236},{-240,-210,-254,-249,-213,-247},{-226,-245,-222,-215,-222,-238},{-238,-249,-238,-224,-213,-245},{-213,-240,-245,-217,-231,-224},{-210,-231,-213,-226,-219,-243},{-210,-245,-231,-236,-245,-236},{-231,-224,-217,-236,-243,-240},{-240,-217,-229,-222,-247,-252},{-247,-254,-231,-213,-238,-238}},{{-222,-217,-247,-229,-238,-215},{-213,-217,-219,-238,-231,-245},{-229,-247,-231,-215,-210,-252},{-233,-219,-236,-222,-233,-254},{-215,-254,-222,-233,-224,-222},{-224,-252,-238,-217,-215,-231},{-233,-226,-210,-243,-233,-219},{-224,-254,-249,-233,-226,-217},{-245,-240,-249,-231,-247,-254},{-249,-240,-233,-222,-236,-226}},{{-245,-247,-247,-226,-245,-231},{-236,-236,-249,-245,-249,-213},{-226,-252,-252,-231,-247,-219},{-224,-233,-249,-226,-217,-252},{-233,-213,-240,-245,-229,-215},{-224,-247,-222,-217,-233,-236},{-236,-245,-236,-229,-208,-217},{-226,-215,-229,-245,-252,-247},{-229,-243,-215,-219,-245,-217},{-238,-252,-233,-229,-224,-219}}},{{{-226,-217,-217,-202,-234,-221},{-234,-215,-200,-221,-228,-221},{-215,-228,-204,-221,-217,-239},{-195,-202,-239,-228,-224,-224},{-226,-217,-197,-204,-206,-239},{-219,-200,-213,-234,-195,-224},{-237,-206,-234,-197,-204,-204},{-226,-210,-239,-226,-195,-213},{-202,-224,-226,-210,-210,-215},{-226,-208,-239,-217,-208,-230}},{{-202,-197,-232,-226,-208,-197},{-215,-210,-213,-221,-221,-210},{-197,-221,-230,-200,-224,-221},{-200,-221,-232,-195,-208,-232},{-208,-234,-195,-208,-195,-197},{-239,-206,-221,-206,-224,-226},{-213,-228,-217,-230,-206,-228},{-234,-206,-213,-239,-206,-224},{-230,-208,-195,-234,-200,-210},{-215,-221,-221,-206,-213,-200}},{{-202,-204,-224,-232,-197,-195},{-228,-215,-219,-206,-239,-217},{-224,-224,-232,-195,-217,-237},{-237,-215,-195,-204,-232,-202},{-219,-239,-237,-208,-195,-200},{-230,-221,-221,-221,-221,-206},{-210,-202,-200,-195,-210,-221},{-237,-234,-195,-237,-215,-221},{-197,-224,-195,-197,-200,-237},{-234,-197,-213,-219,-224,-221}}},{{{-279,-293,-247,-271,-242,-261},{-287,-263,-247,-293,-255,-277},{-261,-245,-287,-277,-253,-279},{-250,-239,-261,-261,-274,-277},{-269,-285,-255,-242,-245,-269},{-253,-279,-279,-282,-287,-285},{-285,-293,-279,-245,-261,-290},{-282,-287,-282,-266,-247,-287},{-271,-293,-277,-274,-269,-269},{-250,-263,-282,-250,-253,-263}},{{-247,-239,-258,-255,-239,-271},{-271,-293,-282,-239,-263,-258},{-293,-239,-258,-242,-271,-290},{-287,-293,-282,-279,-245,-290},{-250,-242,-261,-285,-282,-269},{-263,-261,-274,-245,-271,-255},{-263,-239,-282,-255,-277,-261},{-253,-274,-255,-258,-271,-239},{-271,-261,-290,-285,-269,-245},{-274,-250,-247,-247,-279,-261}},{{-287,-269,-290,-266,-285,-245},{-245,-293,-239,-242,-266,-261},{-255,-258,-285,-271,-239,-239},{-266,-279,-250,-250,-282,-255},{-250,-239,-293,-274,-263,-263},{-274,-247,-253,-242,-282,-266},{-239,-282,-266,-245,-287,-261},{-255,-287,-290,-282,-263,-250},{-271,-253,-293,-253,-266,-255},{-261,-287,-242,-253,-274,-258}}},{{{-300,-260,-272,-252,-263,-280},{-255,-274,-294,-283,-274,-258},{-252,-291,-260,-266,-294,-277},{-252,-255,-300,-288,-269,-286},{-297,-308,-286,-283,-302,-288},{-258,-274,-269,-297,-283,-288},{-300,-302,-269,-260,-266,-252},{-291,-263,-308,-308,-277,-286},{-272,-300,-277,-258,-294,-294},{-305,-258,-302,-288,-305,-283}},{{-294,-308,-297,-277,-291,-274},{-277,-258,-302,-308,-300,-291},{-294,-302,-252,-283,-263,-266},{-305,-272,-272,-305,-269,-294},{-255,-263,-255,-297,-272,-297},{-283,-305,-263,-274,-260,-252},{-297,-258,-283,-286,-274,-269},{-272,-291,-266,-305,-277,-260},{-308,-274,-288,-297,-308,-269},{-266,-277,-263,-297,-300,-283}},{{-291,-297,-288,-252,-258,-294},{-288,-258,-252,-255,-288,-294},{-252,-252,-252,-263,-291,-302},{-260,-302,-258,-255,-300,-280},{-255,-302,-266,-269,-286,-305},{-300,-286,-300,-269,-302,-300},{-274,-308,-294,-291,-294,-266},{-260,-305,-286,-286,-252,-283},{-263,-272,-288,-302,-294,-274},{-258,-255,-291,-272,-286,-308}}},{{{-300,-260,-294,-266,-291,-308},{-288,-286,-283,-302,-277,-252},{-272,-258,-252,-294,-300,-288},{-260,-305,-277,-260,-269,-294},{-308,-294,-286,-288,-266,-297},{-260,-286,-280,-274,-260,-272},{-308,-272,-291,-286,-277,-272},{-266,-291,-272,-288,-252,-269},{-274,-294,-255,-288,-302,-286},{-252,-283,-260,-286,-302,-308}},{{-283,-255,-263,-260,-255,-255},{-291,-283,-305,-269,-269,-305},{-283,-288,-308,-266,-252,-288},{-300,-300,-255,-269,-300,-266},{-288,-280,-277,-283,-288,-272},{-280,-302,-269,-269,-286,-297},{-277,-283,-302,-297,-308,-255},{-252,-291,-308,-305,-269,-297},{-294,-277,-288,-288,-255,-291},{-294,-286,-272,-277,-297,-297}},{{-302,-272,-255,-305,-297,-255},{-280,-291,-269,-274,-286,-258},{-291,-274,-280,-272,-252,-269},{-283,-255,-272,-260,-300,-255},{-300,-266,-280,-286,-266,-302},{-294,-272,-258,-260,-288,-277},{-263,-258,-277,-255,-266,-280},{-291,-280,-297,-277,-266,-269},{-260,-283,-288,-269,-305,-258},{-252,-269,-305,-302,-294,-302}}},{{{-291,-266,-277,-263,-258,-260},{-280,-308,-302,-274,-302,-283},{-260,-291,-255,-280,-305,-277},{-297,-288,-269,-263,-286,-277},{-272,-283,-260,-277,-272,-260},{-291,-305,-288,-291,-291,-272},{-272,-255,-302,-280,-283,-302},{-277,-274,-305,-283,-258,-283},{-300,-305,-283,-302,-286,-252},{-269,-255,-283,-252,-297,-291}},{{-260,-294,-283,-283,-286,-305},{-291,-286,-302,-252,-302,-294},{-297,-263,-308,-266,-255,-286},{-274,-297,-260,-305,-300,-269},{-294,-305,-263,-272,-283,-263},{-274,-263,-297,-297,-266,-280},{-274,-305,-266,-291,-302,-252},{-294,-274,-288,-260,-283,-297},{-305,-266,-280,-274,-308,-302},{-305,-266,-286,-300,-263,-302}},{{-288,-266,-294,-283,-294,-308},{-300,-274,-258,-255,-274,-260},{-263,-283,-297,-263,-280,-277},{-308,-291,-263,-272,-260,-280},{-260,-308,-263,-280,-255,-286},{-263,-302,-258,-252,-280,-258},{-263,-308,-294,-277,-308,-308},{-286,-308,-291,-263,-255,-297},{-302,-255,-286,-286,-263,-288},{-294,-300,-263,-308,-260,-297}}} }; // supply (positive for supply and negative for demand) under different scenarios 
int iniv[I] = { 2070,261,23,150,170,304,400,300,160 }; // initial inventory 
int ubiv[I] = { 2740,522,414,300,341,608,800,600,320 }; // the upper bounds of inventory 
int len[I1][I2] = { {3,6,2,2,3,3,2,1} }; // the traveling time between two sites 
// cost parameters: 
double c1[I] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 }; // inventory holding cost 
double c2[I] = { 0.05,7.19,13.32,4.32,2.64,5.16,6.85,4.63,1.82 }; // cost of leftover and unmet 
double c3[BN] = { 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 }; // transportation cost for capacity contracts 
double c4[I1][I2] = { {6.54,12.11,3.93,2.4,4.69,6.23,4.21,1.65} }; // cost of spot market 
// bid parameters: 
int RLBN[I1][I2] = { {9,9,9,9,9,9,9,9} }; // number of bids on each arc 
int LBset[I1][I2][LBN] = { {{0,1,2,3,4,5,6,7,8},{9,10,11,12,13,14,15,16,17},{18,19,20,21,22,23,24,25,26},{27,28,29,30,31,32,33,34,35},{36,37,38,39,40,41,42,43,44},{45,46,47,48,49,50,51,52,53},{54,55,56,57,58,59,60,61,62},{63,64,65,66,67,68,69,70,71}} }; // set of bids on each arc 
double frt[BN] = { 5.23,5.23,5.23,4.58,4.58,4.58,3.92,3.92,3.92,9.69,9.69,9.69,8.48,8.48,8.48,7.27,7.27,7.27,3.14,3.14,3.14,2.75,2.75,2.75,2.36,2.36,2.36,1.92,1.92,1.92,1.68,1.68,1.68,1.44,1.44,1.44,3.75,3.75,3.75,3.28,3.28,3.28,2.81,2.81,2.81,4.98,4.98,4.98,4.36,4.36,4.36,3.74,3.74,3.74,3.37,3.37,3.37,2.95,2.95,2.95,2.53,2.53,2.53,1.32,1.32,1.32,1.16,1.16,1.16,0.99,0.99,0.99 }; // freight rate of each bid 
int lbcap[BN] = { 75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226,75,75,75,151,151,151,226,226,226}; // capacity lower bound of each bid  
int ubcap[BN] = { 150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300,150,150,150,225,225,225,300,300,300 }; // capacity upper bound of each bid 
// shipment parameters: 
int SHN[BN] = { 8,4,3,7,4,3,8,4,3,6,3,2,5,3,2,6,3,2,7,4,3,7,4,3,8,4,3,7,4,3,7,4,3,8,4,3,8,4,3,7,4,3,8,4,3,8,4,3,7,4,3,7,4,3,7,4,3,7,4,3,8,4,3,9,4,3,9,4,3,8,4,3 }; // number of shipment for each bid 
int SHset[BN][MBSN] = { {0,1,2,3,4,5,6,7},{8,9,10,11},{12,13,14},{15,16,17,18,19,20,21},{22,23,24,25},{26,27,28},{29,30,31,32,33,34,35,36},{37,38,39,40},{41,42,43},{44,45,46,47,48,49},{50,51,52},{53,54},{55,56,57,58,59},{60,61,62},{63,64},{65,66,67,68,69,70},{71,72,73},{74,75},{76,77,78,79,80,81,82},{83,84,85,86},{87,88,89},{90,91,92,93,94,95,96},{97,98,99,100},{101,102,103},{104,105,106,107,108,109,110,111},{112,113,114,115},{116,117,118},{119,120,121,122,123,124,125},{126,127,128,129},{130,131,132},{133,134,135,136,137,138,139},{140,141,142,143},{144,145,146},{147,148,149,150,151,152,153,154},{155,156,157,158},{159,160,161},{162,163,164,165,166,167,168,169},{170,171,172,173},{174,175,176},{177,178,179,180,181,182,183},{184,185,186,187},{188,189,190},{191,192,193,194,195,196,197,198},{199,200,201,202},{203,204,205},{206,207,208,209,210,211,212,213},{214,215,216,217},{218,219,220},{221,222,223,224,225,226,227},{228,229,230,231},{232,233,234},{235,236,237,238,239,240,241},{242,243,244,245},{246,247,248},{249,250,251,252,253,254,255},{256,257,258,259},{260,261,262},{263,264,265,266,267,268,269},{270,271,272,273},{274,275,276},{277,278,279,280,281,282,283,284},{285,286,287,288},{289,290,291},{292,293,294,295,296,297,298,299,300},{301,302,303,304},{305,306,307},{308,309,310,311,312,313,314,315,316},{317,318,319,320},{321,322,323},{324,325,326,327,328,329,330,331},{332,333,334,335},{336,337,338} }; // set of shipments in each bid 
int SHsts[SPN] = { 0,2,4,6,8,10,12,14,2,6,10,14,1,7,13,2,4,6,8,10,12,14,2,6,10,14,0,6,12,0,2,4,6,8,10,12,14,2,6,10,14,0,6,12,0,2,4,6,8,10,2,6,10,1,7,2,4,6,8,10,0,4,8,2,8,1,3,5,7,9,11,1,5,9,1,7,2,4,6,8,10,12,14,0,4,8,12,1,7,13,2,4,6,8,10,12,14,1,5,9,13,0,6,12,0,2,4,6,8,10,12,14,2,6,10,14,0,6,12,2,4,6,8,10,12,14,1,5,9,13,1,7,13,2,4,6,8,10,12,14,0,4,8,12,1,7,13,0,2,4,6,8,10,12,14,2,6,10,14,0,6,12,0,2,4,6,8,10,12,14,1,5,9,13,0,6,12,2,4,6,8,10,12,14,2,6,10,14,2,8,14,0,2,4,6,8,10,12,14,0,4,8,12,2,8,14,0,2,4,6,8,10,12,14,0,4,8,12,1,7,13,2,4,6,8,10,12,14,2,6,10,14,1,7,13,2,4,6,8,10,12,14,0,4,8,12,2,8,14,2,4,6,8,10,12,14,0,4,8,12,1,7,13,2,4,6,8,10,12,14,1,5,9,13,0,6,12,1,3,5,7,9,11,13,15,1,5,9,13,1,7,13,0,2,4,6,8,10,12,14,16,1,5,9,13,0,6,12,0,2,4,6,8,10,12,14,16,1,5,9,13,1,7,13,2,4,6,8,10,12,14,16,1,5,9,13,2,8,14 }; // set of start times of each shipment 
int SHets[SPN] = { 3,5,7,9,11,13,15,17,5,9,13,17,4,10,16,5,7,9,11,13,15,17,5,9,13,17,3,9,15,3,5,7,9,11,13,15,17,5,9,13,17,3,9,15,6,8,10,12,14,16,8,12,16,7,13,8,10,12,14,16,6,10,14,8,14,7,9,11,13,15,17,7,11,15,7,13,4,6,8,10,12,14,16,2,6,10,14,3,9,15,4,6,8,10,12,14,16,3,7,11,15,2,8,14,2,4,6,8,10,12,14,16,4,8,12,16,2,8,14,4,6,8,10,12,14,16,3,7,11,15,3,9,15,4,6,8,10,12,14,16,2,6,10,14,3,9,15,2,4,6,8,10,12,14,16,4,8,12,16,2,8,14,3,5,7,9,11,13,15,17,4,8,12,16,3,9,15,5,7,9,11,13,15,17,5,9,13,17,5,11,17,3,5,7,9,11,13,15,17,3,7,11,15,5,11,17,3,5,7,9,11,13,15,17,3,7,11,15,4,10,16,5,7,9,11,13,15,17,5,9,13,17,4,10,16,5,7,9,11,13,15,17,3,7,11,15,5,11,17,4,6,8,10,12,14,16,2,6,10,14,3,9,15,4,6,8,10,12,14,16,3,7,11,15,2,8,14,3,5,7,9,11,13,15,17,3,7,11,15,3,9,15,1,3,5,7,9,11,13,15,17,2,6,10,14,1,7,13,1,3,5,7,9,11,13,15,17,2,6,10,14,2,8,14,3,5,7,9,11,13,15,17,2,6,10,14,3,9,15 }; // set of end times of each shipment 

//input data ends here

int TSN;
vector<double> trpr;

extern const int ND1 = I1 * T;
extern const int ND2 = I2 * T;
extern const int ND = I * T;
int Nodeset[2][ND];
int NDindex[I][T];
double COST1[ND];
double COST2[ND];
vector<vector<double>> D;
double TQ[ND];
double IQ[ND];

vector<vector<int>> ANDset(2);
vector<double> COST3;
int ARC1 = 0; // this is the number of arcs associated with bids
int ARC = 0; // this is the total no. of arcs

int BAN[BN]; // number of arcs associated with each bid
vector<vector<int>> BAset(BN); // set of arcs associated with each bid
int INAN[ND];
vector<vector<int>> INAset(ND);
int OUTAN[ND];
vector<vector<int>> OUTAset(ND);

double mfb[BN];

double Ysol[BN]; // this is the solution for y variables
double BYsol[BN]; // this is the solution for y variables under the best solution
int TBCN = 0; // total number of Benders cuts
vector < vector<double>> BCyval(BN);
vector<double> BCzval;

// the following are for benders
int TCN5; // total number of constraints for 4
int C5index[I1][T]; // index for C5
int TCN7; // total number of constraints for 4
int C7index[I2][T]; // index for C5

double tempbcy[BN];
double tempbcz;
double bidcost;
vector<double> SMD; // the simulated demand
vector<double> MMD; // the minimum demand at each demand node
vector<double> SCD; // the demand at each demand node under each scenario
vector<double> SAV; // the simulated traffic volume
// the following are for sampling
double UB; // upper bound
double LB; // lower bound
double stgap; // statisical gap
double t_total;
double TTspotcost = 0;
double TTivcost = 0;
double TTbgcost = 0;
double TTbidcap = 0;
double TTbidvol = 0;
double TTspotvol = 0;
double TTivvol = 0;
double TTbgvol = 0;
double TTivvol_s = 0;
double TTbgvol_s = 0;
double TTivvol_d = 0;
double TTbgvol_d = 0;

string num2str(double i)
{
	stringstream ss;
	ss << i;
	return ss.str();
}

ofstream mycout01("//scratch//lxwu//2024//BENDERS//Results//D" + num2str(Dev) + "_P" + num2str(P) + "_SN" + num2str(SN) + "_CN" + num2str(CNN) + ".txt");

void Netsetup() {
	int counter1 = 0;
	for (int i = 0; i != I1; ++i) {
		for (int t = 0; t != T; ++t) {
			Nodeset[0][counter1] = i;
			Nodeset[1][counter1] = t;
			COST1[counter1] = c1[i];
			COST2[counter1] = c2[i];
			TQ[counter1] = ubiv[i];
			IQ[counter1] = iniv[i];
			NDindex[i][t] = counter1;
			counter1 += 1;
		}
	}
	for (int i = 0; i != I2; ++i) {
		for (int t = 0; t != T; ++t) {
			Nodeset[0][counter1] = i + I1;
			Nodeset[1][counter1] = t;
			COST1[counter1] = c1[i + I1];
			COST2[counter1] = c2[i + I1];
			TQ[counter1] = ubiv[i + I1];
			IQ[counter1] = iniv[i + I1];
			NDindex[i + I1][t] = counter1;
			counter1 += 1;
		}
	}
	for (int n = 0; n != ND; ++n) {
		INAN[n] = 0;
		OUTAN[n] = 0;
	}
	counter1 = 0;
	for (int i1 = 0; i1 != I1; ++i1) {
		for (int i2 = 0; i2 != I2; ++i2) {
			for (int n = 0; n != RLBN[i1][i2]; ++n) {
				int thebid = LBset[i1][i2][n];
				for (int m = 0; m != SHN[thebid]; ++m) {// number of shipments
					int theship = SHset[thebid][m];
					int thestart = SHsts[theship];  // start period
					int theend = SHets[theship];    //end period
					int nd1 = NDindex[i1][thestart]; // index of the location-time
					int nd2 = NDindex[i2 + I1][theend];
					BAN[thebid] += 1;
					BAset[thebid].resize(BAN[thebid], counter1);
					INAN[nd2] += 1;
					INAset[nd2].resize(INAN[nd2], counter1);
					OUTAN[nd1] += 1;
					OUTAset[nd1].resize(OUTAN[nd1], counter1);
					counter1 += 1;
					ANDset[0].resize(counter1, nd1); // the head of an arc
					ANDset[1].resize(counter1, nd2); // the tail of an arc
					COST3.resize(counter1, c3[thebid]); // the cost of the arc
				}
			}
		}
	}
	ARC1 = counter1; // total number of arcs associated with bids
	for (int i1 = 0; i1 != I1; ++i1) {
		for (int t1 = 0; t1 != T; ++t1) {
			for (int i2 = 0; i2 != I2; ++i2) {
				if (t1 + len[i1][i2] < T) {
					int t2 = t1 + len[i1][i2];
					int nd1 = NDindex[i1][t1];
					int nd2 = NDindex[i2 + I1][t2];
					INAN[nd2] += 1;
					INAset[nd2].resize(INAN[nd2], counter1);
					OUTAN[nd1] += 1;
					OUTAset[nd1].resize(OUTAN[nd1], counter1);
					counter1 += 1;
					ANDset[0].resize(counter1, nd1);
					ANDset[1].resize(counter1, nd2);
					COST3.resize(counter1, c4[i1][i2]);
				}
			}
		}
	}
	ARC = counter1;
}

void Benderssetup() {
	int counter = 0;
	for (int i = 0; i != I1; ++i) {
		for (int t = 1; t != T; ++t) {
			C5index[i][t] = counter; // index for C5
			counter += 1;
		}
	}
	TCN5 = counter; // number of constraints 5
	counter = 0;
	for (int i = 0; i != I2; ++i) {
		for (int t = 1; t != T; ++t) {
			C7index[i][t] = counter; // index for C7
			counter += 1;
		}
	}
	TCN7 = counter; // number of constraints 7
}

void Treesetup() {
	TSN = 0;
	int search[P];
	for (int p = 0; p != P - 1; ++p) {
		search[p] = 0;
	}
	int counter = 0;
	int doable = 1;
	while (doable == 1) {
		for (int s = 0; s != SN; ++s) {
			search[P - 1] = s;
			double temppro = PRO[P - 1][s];
			D.resize(counter + 1);
			D[counter].resize(ND);
			for (int p = 0; p != P - 1; ++p) {
				temppro = temppro * PRO[p][search[p]];
				for (int n = 0; n != PTN[p]; ++n) {
					int thetime = PTset[p][n];
					for (int i = 0; i != I; ++i) {
						int thend = NDindex[i][thetime];
						D[counter][thend] = sup[i][p][search[p]][n];
					}
				}
			}
			for (int n = 0; n != PTN[P - 1]; ++n) {
				int thetime = PTset[P - 1][n];
				for (int i = 0; i != I; ++i) {
					int thend = NDindex[i][thetime];
					D[counter][thend] = sup[i][P - 1][s][n];
				}
			}
			trpr.resize(counter + 1, temppro);
			counter += 1;
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
	TSN = counter; // this is the total number of scenarios
}

void Rdsample(int M) {
	TSN = M;
	for (int nn = 0; nn != TSN; ++nn) {
		D.resize(nn + 1);
		D[nn].resize(ND);
		for (int p = 0; p != P; ++p) {
			int thess = random(SN); // randomly generate the sample scenario
			for (int n = 0; n != PTN[p]; ++n) {
				int thetime = PTset[p][n];
				for (int i = 0; i != I; ++i) {
					int thend = NDindex[i][thetime];
					D[nn][thend] = sup[i][p][thess][n];
				}
			}
		}
		trpr.resize(nn + 1, double(1) / double(M));
	}
}

void Sampling() {
	if (pow(SN, P) <= 1000) {
		Treesetup();
	}
	else {
		Rdsample(1000);
	}
}

void Cplexsetup() {
	for (int b = 0; b != BN; ++b) {
		mfb[b] = frt[b] * SHN[b];
	}
}

void Simusetup() { // to set up the simulation
	SMD.resize(ND, 0);
	MMD.resize(ND, 0);
	SCD.resize(ND, 0);
	for (int i = 0; i != I1; ++i) {
		for (int p = 0; p != P; ++p) {
			for (int tn = 0; tn != PTN[p]; ++tn) {
				int thet = PTset[p][tn];
				int ndindex = NDindex[i][thet];
				double thedem = 0;
				double themd = 9999999;
				for (int s = 0; s != SN; ++s) {
					thedem += PRO[p][s] * sup[i][p][s][tn];
					if (themd > sup[i][p][s][tn]) { // find the secnario with mini demand (demand are in negtive terms at demand locations)
						themd = sup[i][p][s][tn];
					}
				}
				MMD[ndindex] = themd;
				SMD[ndindex] = thedem;
			}
		}
	}
	for (int i = 0; i != I2; ++i) {
		for (int p = 0; p != P; ++p) {
			for (int tn = 0; tn != PTN[p]; ++tn) {
				int thet = PTset[p][tn];
				int ndindex = NDindex[i + I1][thet];
				double thedem = 0;
				double themd = -9999999;
				for (int s = 0; s != SN; ++s) {
					thedem += PRO[p][s] * sup[i + I1][p][s][tn];
					if (themd < sup[i + I1][p][s][tn]) { // find the secnario with mini demand (demand are in negtive terms at demand locations)
						themd = sup[i + I1][p][s][tn];
					}
				}
				MMD[ndindex] = themd;
				SMD[ndindex] = thedem;
			}
		}
	}
	for (int i = 0; i != I; ++i) {
		for (int tn = 0; tn != PTN[0]; ++tn) {
			int thet = PTset[0][tn];
			int ndindex = NDindex[i][thet];
			MMD[ndindex] = SMD[ndindex];
		}
	}
	SAV.resize(ARC, 0);
}

void inputData() {
	Netsetup();
	Sampling();
	Cplexsetup();
	Benderssetup();
	Simusetup();
}

double Master() {
	double thecost = 0;
	ILOSTLBEGIN
	IloEnv env;
	IloCplex cplex(env);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Param::TimeLimit, 9000.0);
	cplex.setParam(IloCplex::Param::Emphasis::Numerical, 1);
	cplex.setParam(IloCplex::Param::WorkMem, 16384.0);
	cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 16384.0);
	cplex.setParam(IloCplex::Threads, 1);
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVarArray x(env, BN, 0, 1, ILOINT);
	IloNumVarArray y(env, BN, 0, IloInfinity, ILOFLOAT);
	IloNumVar eta(env, 0, IloInfinity, ILOFLOAT);
	try {
		IloInt b;
		IloModel model(env);
		{
			model.add(IloMinimize(env, obj));
		}
		{
			IloExpr Z1(env);
			for (b = 0; b != BN; ++b) {
				Z1 += mfb[b] * y[b];
			}
			model.add(obj == Z1 + eta);
			Z1.end();
		}
		{
			for (b = 0; b != BN; ++b) {
				model.add(y[b] >= lbcap[b] * x[b]);
			}
		}
		{
			for (b = 0; b != BN; ++b) {
				model.add(y[b] <= ubcap[b] * x[b]);
			}
		}
		{ // Benders cut:
			for (int cn = 0; cn != TBCN; ++cn) {
				IloExpr sum1(env);
				for (int bn = 0; bn != BN; ++bn) {
					sum1 += y[bn] * BCyval[bn][cn];
				}
				model.add(eta >= sum1 + BCzval[cn]);
				sum1.end();
			}
		}
		{
			cplex.extract(model);
			cplex.solve();
			//cout << "The model gap is" << cplex.getMIPRelativeGap() * 100 << endl;
			thecost = cplex.getObjValue();
			//cout << "The result is" << thecost << endl;
			bidcost = 0;
			for (b = 0; b != BN; ++b) {
				Ysol[b] = 0;
				double theres = cplex.getValue(x[b]);
				if (theres > 0.99) {
					double theamount = cplex.getValue(y[b]);
					bidcost += mfb[b] * theamount;
					Ysol[b] = theamount; // the solution to y
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
	return thecost;
}

double Subproblem(int sn) { // sn is the scenario
	for (int i = 0; i != I; ++i) {
		for (int tn = 0; tn != PTN[0]; ++tn) {
			int thet = PTset[0][tn];
			int ndindex = NDindex[i][thet];
			MMD[ndindex] = D[sn][ndindex];
		}
	}
	double thecost = 0;
	ILOSTLBEGIN
		IloEnv env;
	IloCplex cplex(env);
	cplex.setOut(env.getNullStream());
	IloNumVarArray var4(env, ARC1, -IloInfinity, 0, ILOFLOAT); // variables associated with constraints 4
	IloNumVarArray var5(env, TCN5, -IloInfinity, IloInfinity, ILOFLOAT); // variables associated with constraints 5
	IloNumVarArray svar5(env, TCN5, -IloInfinity, IloInfinity, ILOFLOAT); // variables associated with constraints 5
	IloNumVarArray var6(env, I1, -IloInfinity, IloInfinity, ILOFLOAT); // variables associated with constraints 6
	IloNumVarArray svar6(env, I1, -IloInfinity, IloInfinity, ILOFLOAT); // variables associated with constraints 6
	IloNumVarArray var7(env, TCN7, -IloInfinity, IloInfinity, ILOFLOAT); // variables associated with constraints 7
	IloNumVarArray svar7(env, TCN7, -IloInfinity, IloInfinity, ILOFLOAT); // variables associated with constraints 7
	IloNumVarArray var8(env, I2, -IloInfinity, IloInfinity, ILOFLOAT); // variables associated with constraints 8
	IloNumVarArray svar8(env, I2, -IloInfinity, IloInfinity, ILOFLOAT); // variables associated with constraints 8
	IloNumVarArray var9(env, ND, -IloInfinity, 0, ILOFLOAT); // variables associated with constraints 9
	IloNumVarArray svar9(env, ND, -IloInfinity, 0, ILOFLOAT); // variables associated with constraints 9
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT); // this is the objective function
	try {
		IloInt b;
		IloModel model(env);
		{
			model.add(IloMaximize(env, obj));
		}
		{
			IloExpr Z1(env);
			IloExpr Z2(env); //D[counter][thend]
			IloExpr Z3(env);
			IloExpr Z4(env);
			IloExpr Z5(env);
			IloExpr Z6(env);
			for (b = 0; b != BN; ++b) { //constraint 4
				for (int an = 0; an != BAN[b]; ++an) {
					int thearc = BAset[b][an]; //locate the arc
					Z1 += Ysol[b] * var4[thearc];
				}
			}
			for (int i = 0; i != I1; ++i) {//constraint 5
				for (int t = 1; t != T; ++t) {
					int theind = C5index[i][t];
					int thend = NDindex[i][t];
					Z2 += D[sn][thend] * var5[theind];
					Z2 += MMD[thend] * svar5[theind];
				}
			}
			for (int i = 0; i != I1; ++i) {//constraint 6
				int thend = NDindex[i][0];
				Z3 += (D[sn][thend] + IQ[thend]) * var6[i];
				Z3 += (MMD[thend] + IQ[thend]) * svar6[i];
			}
			for (int i = 0; i != I2; ++i) {//constraint 7
				for (int t = 1; t != T; ++t) {
					int theind = C7index[i][t];
					int thend = NDindex[i + I1][t];
					Z4 += D[sn][thend] * var7[theind];
					Z4 += MMD[thend] * svar7[theind];
				}
			}
			for (int i = 0; i != I2; ++i) {//constraint 8
				int thend = NDindex[i + I1][0];
				Z5 += (D[sn][thend] + IQ[thend]) * var8[i];
				Z5 += (MMD[thend] + IQ[thend]) * svar8[i];
			}
			for (int nd = 0; nd != ND; ++nd) {//constraint 9
				Z6 += TQ[nd] * var9[nd];
				Z6 += TQ[nd] * svar9[nd];
			}
			model.add(obj == Z1 + Z2 + Z3 + Z4 + Z5 + Z6);
			Z1.end();
			Z2.end();
			Z3.end();
			Z4.end();
			Z5.end();
			Z6.end();
		}
		{ // constraints for z associated with bids
			for (int a = 0; a != ARC1; ++a) {
				int theoutnd = ANDset[0][a]; // this is the out nd
				int theouti = Nodeset[0][theoutnd];
				int theoutt = Nodeset[1][theoutnd];
				int theinnd = ANDset[1][a]; // this is the out nd
				int theini = Nodeset[0][theinnd] - I1;
				int theint = Nodeset[1][theinnd];
				IloExpr sum1(env);
				sum1 += var4[a];
				if (theoutt == 0) {
					sum1 += var6[theouti] + svar6[theouti];
				}
				else {
					int theind = C5index[theouti][theoutt];
					sum1 += var5[theind] + svar5[theind];
				}
				if (theint == 0) {
					sum1 -= (var8[theini] + svar8[theini]);
				}
				else {
					int theind = C7index[theini][theint];
					sum1 -= (var7[theind] + svar7[theind]);
				}
				model.add(sum1 <= COST3[a]);
				sum1.end();
			}
		}
		{ // constraints for z not associated with bids
			for (int a = ARC1; a != ARC; ++a) {
				int theoutnd = ANDset[0][a]; // this is the out nd
				int theouti = Nodeset[0][theoutnd];
				int theoutt = Nodeset[1][theoutnd];
				int theinnd = ANDset[1][a]; // this is the out nd
				int theini = Nodeset[0][theinnd] - I1;
				int theint = Nodeset[1][theinnd];
				IloExpr sum1(env);
				if (theoutt == 0) {
					sum1 += var6[theouti] + svar6[theouti];
				}
				else {
					int theind = C5index[theouti][theoutt];
					sum1 += var5[theind] + svar5[theind];
				}
				if (theint == 0) {
					sum1 -= (var8[theini] + svar8[theini]);
				}
				else {
					int theind = C7index[theini][theint];
					sum1 -= (var7[theind] + svar7[theind]);
				}
				model.add(sum1 <= COST3[a]);
				sum1.end();
			}
		}
		{ // constraints for u associated with supply nodes
			for (int i = 0; i != I1; ++i) {
				IloExpr sum1(env); // for t=0
				int theind = C5index[i][1];
				sum1 += -var5[theind];
				sum1 += var6[i];
				int thend = NDindex[i][0];
				sum1 += var9[thend];
				model.add(sum1 <= COST1[thend]);
				sum1.end();
				for (int t = 1; t != T; ++t) {
					IloExpr sum1(env);
					int theind1 = C5index[i][t];
					sum1 += var5[theind1];
					if (t < T - 1) {
						int theind2 = C5index[i][t + 1];
						sum1 += -var5[theind2];
					}
					int thend = NDindex[i][t];
					sum1 += var9[thend];
					model.add(sum1 <= COST1[thend]);
					sum1.end();
				}
			}
		}
		{ // constraints for su associated with supply nodes
			for (int i = 0; i != I1; ++i) {
				IloExpr sum1(env); // for t=0
				int theind = C5index[i][1];
				sum1 += -svar5[theind];
				sum1 += svar6[i];
				int thend = NDindex[i][0];
				sum1 += svar9[thend];
				model.add(sum1 <= 0);
				sum1.end();
				for (int t = 1; t != T; ++t) {
					IloExpr sum1(env);
					int theind1 = C5index[i][t];
					sum1 += svar5[theind1];
					if (t < T - 1) {
						int theind2 = C5index[i][t + 1];
						sum1 += -svar5[theind2];
					}
					int thend = NDindex[i][t];
					sum1 += svar9[thend];
					model.add(sum1 <= 0);
					sum1.end();
				}
			}
		}
		{ // constraints for u associated with demand nodes
			for (int i = 0; i != I2; ++i) {
				IloExpr sum1(env); // for t=0
				int theind = C7index[i][1];
				sum1 += -var7[theind];
				sum1 += var8[i];
				int thend = NDindex[i + I1][0];
				sum1 += var9[thend];
				model.add(sum1 <= COST1[thend]);
				sum1.end();
				for (int t = 1; t != T; ++t) {
					IloExpr sum1(env);
					int theind1 = C7index[i][t];
					sum1 += var7[theind1];
					if (t < T - 1) {
						int theind2 = C7index[i][t + 1];
						sum1 += -var7[theind2];
					}
					int thend = NDindex[i + I1][t];
					sum1 += var9[thend];
					model.add(sum1 <= COST1[thend]);
					sum1.end();
				}
			}
		}
		{ // constraints for su associated with demand nodes
			for (int i = 0; i != I2; ++i) {
				IloExpr sum1(env); // for t=0
				int theind = C7index[i][1];
				sum1 += -svar7[theind];
				sum1 += svar8[i];
				int thend = NDindex[i + I1][0];
				sum1 += svar9[thend];
				model.add(sum1 <= 0);
				sum1.end();
				for (int t = 1; t != T; ++t) {
					IloExpr sum1(env);
					int theind1 = C7index[i][t];
					sum1 += svar7[theind1];
					if (t < T - 1) {
						int theind2 = C7index[i][t + 1];
						sum1 += -svar7[theind2];
					}
					int thend = NDindex[i + I1][t];
					sum1 += svar9[thend];
					model.add(sum1 <= 0);
					sum1.end();
				}
			}
		}
		{ // constraints for v associated with supply nodes
			for (int i = 0; i != I1; ++i) {
				IloExpr sum1(env); // for t=0
				int theind = C5index[i][1];
				sum1 += -var5[theind];
				sum1 += var6[i];
				int thend = NDindex[i][0];
				model.add(sum1 <= COST2[thend]);
				sum1.end();
				for (int t = 1; t != T; ++t) {
					IloExpr sum1(env);
					int theind1 = C5index[i][t];
					sum1 += var5[theind1];
					if (t < T - 1) {
						int theind2 = C5index[i][t + 1];
						sum1 += -var5[theind2];
					}
					int thend = NDindex[i][t];
					model.add(sum1 <= COST2[thend]);
					sum1.end();
				}
			}
		}
		{ // constraints for sv associated with supply nodes
			for (int i = 0; i != I1; ++i) {
				IloExpr sum1(env); // for t=0
				int theind = C5index[i][1];
				sum1 += -svar5[theind];
				sum1 += svar6[i];
				int thend = NDindex[i][0];
				model.add(sum1 <= 0);
				sum1.end();
				for (int t = 1; t != T; ++t) {
					IloExpr sum1(env);
					int theind1 = C5index[i][t];
					sum1 += svar5[theind1];
					if (t < T - 1) {
						int theind2 = C5index[i][t + 1];
						sum1 += -svar5[theind2];
					}
					int thend = NDindex[i][t];
					model.add(sum1 <= 0);
					sum1.end();
				}
			}
		}
		{ // constraints for v associated with demand nodes
			for (int i = 0; i != I2; ++i) {
				IloExpr sum1(env); // for t=0
				int theind = C7index[i][1];
				sum1 += var7[theind];
				sum1 += -var8[i];
				int thend = NDindex[i + I1][0];
				model.add(sum1 <= COST2[thend]);
				sum1.end();
				for (int t = 1; t != T; ++t) {
					IloExpr sum1(env);
					int theind1 = C7index[i][t];
					sum1 += -var7[theind1];
					if (t < T - 1) {
						int theind2 = C7index[i][t + 1];
						sum1 += var7[theind2];
					}
					int thend = NDindex[i + I1][t];
					model.add(sum1 <= COST2[thend]);
					sum1.end();
				}
			}
		}
		{ // constraints for sv associated with demand nodes
			for (int i = 0; i != I2; ++i) {
				IloExpr sum1(env); // for t=0
				int theind = C7index[i][1];
				sum1 += svar7[theind];
				sum1 += -svar8[i];
				int thend = NDindex[i + I1][0];
				model.add(sum1 <= 0);
				sum1.end();
				for (int t = 1; t != T; ++t) {
					IloExpr sum1(env);
					int theind1 = C7index[i][t];
					sum1 += -svar7[theind1];
					if (t < T - 1) {
						int theind2 = C7index[i][t + 1];
						sum1 += svar7[theind2];
					}
					int thend = NDindex[i + I1][t];
					model.add(sum1 <= 0);
					sum1.end();
				}
			}
		}
		{
			cplex.extract(model);
			cplex.solve();
			for (b = 0; b != BN; ++b) { //constraint 4
				double theytemp = 0;
				for (int an = 0; an != BAN[b]; ++an) {
					int thearc = BAset[b][an]; //locate the arc
					theytemp += cplex.getValue(var4[thearc]);
				}
				tempbcy[b] += trpr[sn] * theytemp;
			}
			double theztemp = 0;
			for (int i = 0; i != I1; ++i) {//constraint 5
				for (int t = 1; t != T; ++t) {
					int theind = C5index[i][t];
					int thend = NDindex[i][t];
					theztemp += D[sn][thend] * cplex.getValue(var5[theind]);
					theztemp += MMD[thend] * cplex.getValue(svar5[theind]);
				}
			}
			for (int i = 0; i != I1; ++i) {//constraint 6
				int thend = NDindex[i][0];
				theztemp += (D[sn][thend] + IQ[thend]) * cplex.getValue(var6[i]);
				theztemp += (MMD[thend] + IQ[thend]) * cplex.getValue(svar6[i]);
			}
			for (int i = 0; i != I2; ++i) {//constraint 7
				for (int t = 1; t != T; ++t) {
					int theind = C7index[i][t];
					int thend = NDindex[i + I1][t];
					theztemp += D[sn][thend] * cplex.getValue(var7[theind]);
					theztemp += MMD[thend] * cplex.getValue(svar7[theind]);
				}
			}
			for (int i = 0; i != I2; ++i) {//constraint 8
				int thend = NDindex[i + I1][0];
				theztemp += (D[sn][thend] + IQ[thend]) * cplex.getValue(var8[i]);
				theztemp += (MMD[thend] + IQ[thend]) * cplex.getValue(svar8[i]);
			}
			for (int nd = 0; nd != ND; ++nd) {//constraint 9
				theztemp += TQ[nd] * cplex.getValue(var9[nd]);
				theztemp += TQ[nd] * cplex.getValue(svar9[nd]);
			}
			tempbcz += trpr[sn] * theztemp;
			thecost = cplex.getObjValue();
		}
	}
	catch (IloException& ex) {
		cerr << ex << endl;
	}
	catch (...) {
		cerr << "Error..." << endl;
	}
	env.end();
	//cout << "The dual result is" << thecost << endl;
	return thecost;
}

void Benders() {
	double thegap = 100;
	double upperbound = 99999999;
	struct timeval t_start, t_crt;   // times
	gettimeofday(&t_start, NULL);
	double timelimt = 1800 * P;
	while (thegap > 0.01) {
		gettimeofday(&t_crt, NULL);
		double thedure = ((t_crt.tv_sec - t_start.tv_sec) * 1000000 + (t_crt.tv_usec - t_start.tv_usec)) / double(1000000);
		if (thedure >= timelimt) {
			break;
		}
		double fcost = Master(); // first stage cost
		double scost = bidcost;
		for (int bn = 0; bn != BN; ++bn) {
			tempbcy[bn] = 0;
		}
		tempbcz = 0;
		for (int sn = 0; sn != TSN; ++sn) {
			for (int i = 0; i != I; ++i) {
				for (int tn = 0; tn != PTN[0]; ++tn) {
					int thet = PTset[0][tn];
					int ndindex = NDindex[i][thet];
					MMD[ndindex] = D[sn][ndindex]; // this is the demand we know at tge second stage
				}
			}
			double subcost = Subproblem(sn);
			scost += trpr[sn] * subcost;
		}
		if (upperbound > scost) {
			upperbound = scost;
			for (int bn = 0; bn != BN; ++bn) {
				BYsol[bn] = 0;
				if (Ysol[bn] > 0.1) {
					BYsol[bn] = Ysol[bn];
				}
			}
		}
		thegap = ((upperbound - fcost) / upperbound) * 100;
		cout << "Upper bound=" << scost << endl;
		cout << "Lower bound=" << fcost << endl;
		cout << "number of Benders Cuts=" << TBCN << endl;
		cout << "thegap=" << thegap << endl;
		TBCN += 1;
		for (int bn = 0; bn != BN; ++bn) {
			BCyval[bn].resize(TBCN, tempbcy[bn]);
		}
		BCzval.resize(TBCN, tempbcz);
		LB = fcost; // update lower bound
	}
}

void Primal() { // feasibility constraints will be added 
	ILOSTLBEGIN
		for (int i = 0; i != I; ++i) {
			for (int tn = 0; tn != PTN[0]; ++tn) {
				int thet = PTset[0][tn];
				int ndindex = NDindex[i][thet];
				SMD[ndindex] = SCD[ndindex]; // this is the demand we know
				MMD[ndindex] = SCD[ndindex]; // this is the demand we know at tge second stage
			}
		}
	IloEnv env;
	IloCplex cplex(env);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVarArray z(env, ARC, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray u(env, ND, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray v(env, ND, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray mu(env, ND, 0, IloInfinity, ILOFLOAT); // this is the inventory at the minimum demand to let them kept under control
	IloNumVarArray mv(env, ND, 0, IloInfinity, ILOFLOAT); // this is the backlog at the minimum demand to let them kept under control
	try {
		IloModel model(env);
		{
			model.add(IloMinimize(env, obj));
		}
		{
			IloExpr Z2(env);
			IloExpr Z3(env);
			IloExpr Z4(env);
			for (int n = 0; n != ND; ++n) {
				Z2 += COST1[n] * u[n];
			}
			for (int n = 0; n != ND; ++n) {
				Z3 += COST2[n] * v[n];
			}
			for (int a = 0; a != ARC; ++a) {
				Z4 += COST3[a] * z[a];
			}
			model.add(obj == Z2 + Z3 + Z4);
			Z2.end();
			Z3.end();
			Z4.end();
		}
		{
			for (int b = 0; b != BN; ++b) {
				for (int n = 0; n != BAN[b]; ++n) {
					int thearc = BAset[b][n];
					model.add(z[thearc] <= BYsol[b]);
				}
			}
		}
		{
			for (int t1 = 1; t1 != T; ++t1) {
				int t2 = t1 - 1;
				for (int i = 0; i != I1; ++i) {
					int nd1 = NDindex[i][t1];
					int nd2 = NDindex[i][t2];
					IloExpr sum1(env);
					for (int n = 0; n != OUTAN[nd1]; ++n) {
						int thearc = OUTAset[nd1][n];
						sum1 += z[thearc];
					}
					model.add(u[nd1] == u[nd2] + SMD[nd1] + v[nd2] - v[nd1] - sum1);
					sum1.end();
				}
			}
		}
		{
			for (int t1 = 1; t1 != T; ++t1) {
				int t2 = t1 - 1;
				for (int i = 0; i != I1; ++i) {
					int nd1 = NDindex[i][t1];
					int nd2 = NDindex[i][t2];
					IloExpr sum1(env);
					for (int n = 0; n != OUTAN[nd1]; ++n) {
						int thearc = OUTAset[nd1][n];
						sum1 += z[thearc];
					}
					model.add(mu[nd1] == mu[nd2] + MMD[nd1] + mv[nd2] - mv[nd1] - sum1);
					sum1.end();
				}
			}
		}
		{
			for (int i = 0; i != I1; ++i) {
				int nd = NDindex[i][0];
				IloExpr sum1(env);
				for (int n = 0; n != OUTAN[nd]; ++n) {
					int thearc = OUTAset[nd][n];
					sum1 += z[thearc];
				}
				model.add(u[nd] == IQ[nd] + SMD[nd] - v[nd] - sum1);
				sum1.end();
			}
		}
		{
			for (int i = 0; i != I1; ++i) {
				int nd = NDindex[i][0];
				IloExpr sum1(env);
				for (int n = 0; n != OUTAN[nd]; ++n) {
					int thearc = OUTAset[nd][n];
					sum1 += z[thearc];
				}
				model.add(mu[nd] == IQ[nd] + MMD[nd] - mv[nd] - sum1);
				sum1.end();
			}
		}
		{
			for (int t1 = 1; t1 != T; ++t1) {
				int t2 = t1 - 1;
				for (int i = 0; i != I2; ++i) {
					int nd1 = NDindex[i + I1][t1];
					int nd2 = NDindex[i + I1][t2];
					IloExpr sum1(env);
					for (int n = 0; n != INAN[nd1]; ++n) {
						int thearc = INAset[nd1][n];
						sum1 += z[thearc];
					}
					model.add(u[nd1] == u[nd2] + SMD[nd1] - v[nd2] + v[nd1] + sum1);
					sum1.end();
				}
			}
		}
		{
			for (int t1 = 1; t1 != T; ++t1) {
				int t2 = t1 - 1;
				for (int i = 0; i != I2; ++i) {
					int nd1 = NDindex[i + I1][t1];
					int nd2 = NDindex[i + I1][t2];
					IloExpr sum1(env);
					for (int n = 0; n != INAN[nd1]; ++n) {
						int thearc = INAset[nd1][n];
						sum1 += z[thearc];
					}
					model.add(mu[nd1] == mu[nd2] + MMD[nd1] - mv[nd2] + mv[nd1] + sum1);
					sum1.end();
				}
			}
		}
		{
			for (int i = 0; i != I2; ++i) {
				int nd = NDindex[i + I1][0];
				IloExpr sum1(env);
				for (int n = 0; n != INAN[nd]; ++n) {
					int thearc = INAset[nd][n];
					sum1 += z[thearc];
				}
				model.add(u[nd] == IQ[nd] + SMD[nd] + v[nd] + sum1);
				sum1.end();
			}
		}
		{
			for (int i = 0; i != I2; ++i) {
				int nd = NDindex[i + I1][0];
				IloExpr sum1(env);
				for (int n = 0; n != INAN[nd]; ++n) {
					int thearc = INAset[nd][n];
					sum1 += z[thearc];
				}
				model.add(mu[nd] == IQ[nd] + MMD[nd] + mv[nd] + sum1);
				sum1.end();
			}
		}
		{
			for (int n = 0; n != ND; ++n) {
				model.add(u[n] <= TQ[n]);
			}
		}
		{
			for (int n = 0; n != ND; ++n) {
				model.add(mu[n] <= TQ[n]);
			}
		}
		{
			cplex.extract(model);
			cplex.solve();
			for (int an = 0; an != ARC; ++an) {
				SAV[an] = cplex.getValue(z[an]);
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
}

double Sceprimal(int sn) { // to solve the scenario based model
	Primal(); //solve the primal problem
	ILOSTLBEGIN
		double thecost = 0;
	IloEnv env;
	IloCplex cplex(env);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	IloNumVar obj(env, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVarArray z(env, ARC, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray u(env, ND, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray v(env, ND, 0, IloInfinity, ILOFLOAT);
	IloNumVarArray mu(env, ND, 0, IloInfinity, ILOFLOAT); // this is the inventory at the minimum demand to let them kept under control
	IloNumVarArray mv(env, ND, 0, IloInfinity, ILOFLOAT); // this is the backlog at the minimum demand to let them kept under control
	try {
		IloModel model(env);
		{
			model.add(IloMinimize(env, obj));
		}
		{
			IloExpr Z2(env);
			IloExpr Z3(env);
			IloExpr Z4(env);
			for (int n = 0; n != ND; ++n) {
				Z2 += COST1[n] * u[n];
			}
			for (int n = 0; n != ND; ++n) {
				Z3 += COST2[n] * v[n];
			}
			for (int a = 0; a != ARC; ++a) {
				Z4 += COST3[a] * z[a];
			}
			model.add(obj == Z2 + Z3 + Z4);
			Z2.end();
			Z3.end();
			Z4.end();
		}
		{
			for (int b = 0; b != BN; ++b) {
				for (int n = 0; n != BAN[b]; ++n) {
					int thearc = BAset[b][n];
					model.add(z[thearc] <= BYsol[b]);
				}
			}
		}
		{
			for (int t1 = 1; t1 != T; ++t1) {
				int t2 = t1 - 1;
				for (int i = 0; i != I1; ++i) {
					int nd1 = NDindex[i][t1];
					int nd2 = NDindex[i][t2];
					IloExpr sum1(env);
					for (int n = 0; n != OUTAN[nd1]; ++n) {
						int thearc = OUTAset[nd1][n];
						sum1 += z[thearc];
					}
					model.add(u[nd1] == u[nd2] + SCD[nd1] + v[nd2] - v[nd1] - sum1);
					sum1.end();
				}
			}
		}
		{
			for (int i = 0; i != I1; ++i) {
				int nd = NDindex[i][0];
				IloExpr sum1(env);
				for (int n = 0; n != OUTAN[nd]; ++n) {
					int thearc = OUTAset[nd][n];
					sum1 += z[thearc];
				}
				model.add(u[nd] == IQ[nd] + SCD[nd] - v[nd] - sum1);
				sum1.end();
			}
		}
		{
			for (int t1 = 1; t1 != T; ++t1) {
				int t2 = t1 - 1;
				for (int i = 0; i != I2; ++i) {
					int nd1 = NDindex[i + I1][t1];
					int nd2 = NDindex[i + I1][t2];
					IloExpr sum1(env);
					for (int n = 0; n != INAN[nd1]; ++n) {
						int thearc = INAset[nd1][n];
						sum1 += z[thearc];
					}
					model.add(u[nd1] == u[nd2] + SCD[nd1] - v[nd2] + v[nd1] + sum1);
					sum1.end();
				}
			}
		}
		{
			for (int i = 0; i != I2; ++i) {
				int nd = NDindex[i + I1][0];
				IloExpr sum1(env);
				for (int n = 0; n != INAN[nd]; ++n) {
					int thearc = INAset[nd][n];
					sum1 += z[thearc];
				}
				model.add(u[nd] == IQ[nd] + SCD[nd] + v[nd] + sum1);
				sum1.end();
			}
		}
		{
			for (int n = 0; n != ND; ++n) {
				model.add(u[n] <= TQ[n]);
			}
		}
		{
			for (int an = 0; an != ARC; ++an) {
				model.add(z[an] == SAV[an]); // fix the transport volume
			}
		}
		{
			cplex.extract(model);
			cplex.solve();
			thecost = cplex.getObjValue();
			for (int nn = 0; nn != ND; ++nn) {
				double theval = cplex.getValue(u[nn]);
				TTivvol += theval;
				if (Nodeset[0][nn] < I1) {
					TTivvol_s += theval;
				}
				else {
					TTivvol_d += theval;
				}
				TTivcost += theval * COST1[nn];
			}
			for (int nn = 0; nn != ND; ++nn) {
				int thend = nn;
				double theval = cplex.getValue(v[nn]);
				TTbgvol += theval;
				if (Nodeset[0][thend] < I1) {
					TTbgvol_s += theval;
				}
				else {
					TTbgvol_d += theval;
				}
				TTbgcost += theval * COST2[thend];
			}
			for (int nn = 0; nn != ARC; ++nn) {
				int thearc = nn;
				double theval = cplex.getValue(z[nn]);
				if (COST3[thearc] < 0.001) {
					TTbidvol += theval;
				}
				if (COST3[thearc] > 0.001) {
					TTspotvol += theval;
					TTspotcost += theval * COST3[thearc];
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
	return thecost;
}

double SCECAl(int SIND) {
	for (int p = 0; p != P; ++p) {
		int thess = random(SN); // randomly generate the sample scenario
		for (int n = 0; n != PTN[p]; ++n) {
			int thetime = PTset[p][n];
			for (int i = 0; i != I; ++i) {
				int thend = NDindex[i][thetime];
				SCD[thend] = sup[i][p][thess][n];
			}
		}
	}
	double theupper = bidcost;
	double transcost = Sceprimal(SIND);
	theupper += transcost;
	return theupper;
}

void Upperbound(int M) { // sampling m scenarios
	vector<double> Fres(M);
	UB = 0;
	int STSN = 0;
	for (int mm = 0; mm != M; ++mm) {
		double theres = SCECAl(mm);
		UB += theres;
		Fres[STSN] = theres;
		STSN += 1;
	}
	double aveub = double(UB) / double(M);
	double totaldev = 0;
	for (int nn = 0; nn != M; ++nn) {
		totaldev += (double(Fres[nn] - aveub)) * (double(Fres[nn] - aveub));
	}
	double avedev = (totaldev) / (M - 1);
	UB = aveub + 1.96 * pow(avedev, 0.5) / pow(double(M), 0.5);
	stgap = double(UB - LB) / double(UB);
	stgap *= 100;
	for (int b = 0; b != BN; ++b) {
		TTbidcap += BYsol[b] * SHN[b];
	}
	TTspotcost = TTspotcost / double(M);
	TTivcost = TTivcost / double(M);
	TTbgcost = TTbgcost / double(M);
	TTbidvol = TTbidvol / double(M);
	TTspotvol = TTspotvol / double(M);
	TTivvol = TTivvol / double(M);
	TTbgvol = TTbgvol / double(M);
	TTivvol_s = TTivvol_s / double(M);
	TTivvol_d = TTivvol_d / double(M);
	TTbgvol_s = TTbgvol_s / double(M);
	TTbgvol_d = TTbgvol_d / double(M);
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
		for (int s = 0; s != SN; ++s) {
			search[P - 1] = s; // load the last stage scenario
			for (int p = 0; p != P - 1; ++p) { // this circulation does not include the last stage 
				for (int n = 0; n != PTN[p]; ++n) {  // record the supply/demand
					int thetime = PTset[p][n];
					for (int i = 0; i != I; ++i) {
						int thend = NDindex[i][thetime];// thend
						SCD[thend] = sup[i][p][search[p]][n];  // we solve scenario by scenario, each time only one scenario is loaded into the tree (as the first scenario)
					}
				}
			}
			for (int n = 0; n != PTN[P - 1]; ++n) {   // record the supply/demand of last stage
				int thetime = PTset[P - 1][n];
				for (int i = 0; i != I; ++i) {
					int thend = NDindex[i][thetime];// thend
					SCD[thend] = sup[i][P - 1][s][n];
				}
			}
			STSN += 1;
			double theupper = bidcost;
			double scetran = 0;
			scetran = Sceprimal(s); // solve the scenario problem
			theupper += scetran;
			UB += theupper;
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
	stgap = double(aveub - LB) / double(aveub);
	stgap *= 100;
	for (int b = 0; b != BN; ++b) {
		TTbidcap += BYsol[b] * SHN[b];
	}
	TTspotcost = TTspotcost / double(STSN);
	TTivcost = TTivcost / double(STSN);
	TTbgcost = TTbgcost / double(STSN);
	TTbidvol = TTbidvol / double(STSN);
	TTspotvol = TTspotvol / double(STSN);
	TTivvol = TTivvol / double(STSN);
	TTbgvol = TTbgvol / double(STSN);
	TTivvol_s = TTivvol_s / double(STSN);
	TTivvol_d = TTivvol_d / double(STSN);
	TTbgvol_s = TTbgvol_s / double(STSN);
	TTbgvol_d = TTbgvol_d / double(STSN);
}

void Simulation() {
	bidcost = 0;
	for (int bn = 0; bn != BN; ++bn) {
		bidcost += BYsol[bn] * mfb[bn];
	}
	if (pow(SN, P) <= 10000) {
		Enumeration();
	}
	else {
		Upperbound(10000);
	}
}

void documentation() {
	cout << "succeeded" << endl;
	if (pow(SN, P) <= 10000) {
		cout << "Real UB" << endl;
	}
	if (pow(SN, P) > 10000) {
		cout << "Statistical UB" << endl;
	}
	cout << "upperbound=" << UB << endl;
	cout << "totalgap=" << stgap << endl;
	cout << "bidcost=" << bidcost << endl;
	for (int bb = 0; bb != BN; ++bb) {
		if (BYsol[bb] > 1) {
			cout << "bid[" << bb << "] purchased at" << BYsol[bb] << endl;
		}
	}
	cout << "Tbidcap=" << TTbidcap << endl;
	cout << "Tbidcost=" << bidcost << endl;
	cout << "Tspotcost=" << TTspotcost << endl;
	cout << "Tivcost=" << TTivcost << endl;
	cout << "Tbgcost=" << TTbgcost << endl;
	cout << "Tbidvol=" << TTbidvol << endl;
	cout << "Tspotvol=" << TTspotvol << endl;
	cout << "Tivvol=" << TTivvol << endl;
	cout << "Tbgvol=" << TTbgvol << endl;
	cout << "Tivvol_s=" << TTivvol_s << endl;
	cout << "Tivvol_d=" << TTivvol_d << endl;
	cout << "Tbgvol_s=" << TTbgvol_s << endl;
	cout << "Tbgvol_d=" << TTbgvol_d << endl;
	if (pow(SN, P) <= 10000) {
		mycout01 << "Real UB" << endl;
	}
	if (pow(SN, P) > 10000) {
		mycout01 << "Statistical UB" << endl;
	}
	mycout01 << "upperbound=" << endl << UB << endl;
	mycout01 << "totalgap=" << endl << stgap << endl;
	mycout01 << "lb_total=" << endl << LB << endl;
	mycout01 << "t_total=" << endl << t_total << endl;
	mycout01 << "bidcost=" << bidcost << endl;
	for (int bb = 0; bb != BN; ++bb) {
		if (BYsol[bb] > 1) {
			mycout01 << "bid[" << bb << "] purchased at" << BYsol[bb] << endl;
		}
	}
	mycout01 << "Tbidcap=" << TTbidcap << endl;
	mycout01 << "Tbidcost=" << bidcost << endl;
	mycout01 << "Tspotcost=" << TTspotcost << endl;
	mycout01 << "Tivcost=" << TTivcost << endl;
	mycout01 << "Tbgcost=" << TTbgcost << endl;
	mycout01 << "Tbidvol=" << TTbidvol << endl;
	mycout01 << "Tspotvol=" << TTspotvol << endl;
	mycout01 << "Tivvol=" << TTivvol << endl;
	mycout01 << "Tbgvol=" << TTbgvol << endl;
	mycout01 << "Tivvol_s=" << TTivvol_s << endl;
	mycout01 << "Tivvol_d=" << TTivvol_d << endl;
	mycout01 << "Tbgvol_s=" << TTbgvol_s << endl;
	mycout01 << "Tbgvol_d=" << TTbgvol_d << endl;
}

int main(int argc, char** argv) {
	inputData();
	cout << "formulation" << endl;
	struct timeval t_start, t_end; 
	gettimeofday(&t_start,NULL);
	Benders();
	gettimeofday(&t_end,NULL);
	cout << "time=" <<((t_end.tv_sec - t_start.tv_sec)*1000000+(t_end.tv_usec - t_start.tv_usec))/double(1000000)<<endl;
	t_total = ((t_end.tv_sec - t_start.tv_sec) * 1000000 + (t_end.tv_usec - t_start.tv_usec)) / double(1000000);
	Simulation();
	documentation();
	cout << "FULL-CPLEX-succeeded" << endl;
	return 0;
}