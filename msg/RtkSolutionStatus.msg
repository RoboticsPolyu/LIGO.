# RTK ambiguity status associated with the simultaneously published odometry.
std_msgs/Header header

uint8 NO_AMBIGUITY=0
uint8 FLOAT=1
uint8 FIXED=2

uint8 solution_type
string solution_name
string detail_status
uint32 active_float
uint32 active_fixed
bool integer_solution_available
float64 lambda_ratio
