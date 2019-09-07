#Packages
from xlrd import open_workbook
import random
from scipy.stats import entropy

class Problem:

    ''' Parameters:
    m : Number of scenarios (#255)
    n :  Number of Tests (#78)
    U_all : Set of all tests {0, 1, ..., n-1} (remain unchanged)
    U : Set of tests that have not been performed (will change during the algorithm)
    H_all : Set of all scenarios {0, 1, ..., m-1} (remain unchanged)
    H : Set of compatible scenarios (will change during the algorithm)
    outcome: Outcome table, a[s][e] corresponds to feedback of scenario s for test e
    pi : Original probabilities for scenarios (remain unchanged)
    p : Current probabilities for scenarios (will change during the algorithm)
    h : Number of star outcomes for each scnenario
    r : Number of star outcomes for each test
    Te_positive: List of all scenarios with positive feedback for each test
    Te_negative: List of all scenarios with negative feedback for each test
    num_copies_all: Number of copies of scnerios based on their star tests
    num_copies: Number of compatible copies of scenarios (will change during the algorithm)
    '''
    
    def __init__(self, input_sheet):  
        self.m = input_sheet.nrows - 2
        self.n = input_sheet.ncols - 5
        self.U_all = set(range(self.n))
        self.H_all = set(range(self.m))
        self.U = set(range(self.n))
        self.H = set(range(self.m))
        self.outcome = [[]]
        self.pi = [0.0] * self.m
        self.p = [0.0] * self.m
        self.h = [0.0] * self.m
        self.r = [0.0] * self.n
        self.Te_positive = [set() for i in range(self.n)]
        self.Te_negative = [set() for i in range(self.n)]
        self.num_copies_all = [1]*self.m #Number of copies of scnerios based on their star tests
        self.num_copies = [1]*self.m #Number of compatible copies of scenarios
        for s in self.H_all:
            self.outcome.insert(s,[])
            for e in self.U_all:
                self.outcome[s].insert(e, input_sheet.cell(s + 2, e + 2).value) #First 2 columns and 2 rows are not outcomes
    
    #Initializing the probabilities
    def init_probabilities(self, input_sheet, alpha):
        alpha2index = {1: 0, 1.5: 2, 2: 1}
        for s in self.H_all:
            self.pi[s] = input_sheet.cell(s + 2, self.n + 2 + alpha2index[alpha]).value
        print("Lower Bound: ", entropy(self.pi,base=2))
        
    def print_h_and_r(self): #In this function we compute paramaters r and h, as well as their average and maximum
        for s in self.H_all:
            self.h[s] = 0.0
        for e in self.U_all:
            self.r[e] = 0.0
        for s in self.H_all:
            for e in self.U_all:
                if(self.outcome[s][e] == 'u'):
                    self.h[s] = self.h[s] + 1
                    self.r[e] = self.r[e] + 1
                    self.num_copies_all[s] = self.num_copies_all[s] * 2
        print("max_h = ", max(self.h), "average_h = ", sum(self.h)/self.m)
        print("max_r = ", max(self.r), "average_r = ", sum(self.r)/self.n)
        
    #Computing T_e^+ and T_e^- for all tests
    def compute_Te_positive(self): #This function compute T_e^+ for all tests e 
        for e in self.U_all:
            self.Te_positive[e] = set([s for s in self.H_all if self.outcome[s][e] == 1])
    
    def compute_Te_negative(self):  #This function compute T_e^- for all tests e
        for e in self.U_all:
            self.Te_negative[e] = set([s for s in self.H_all if self.outcome[s][e] == 0])

#Main body of Algorithm
def one_step_ODT(this_problem): #This function computes the score of each test and find the test with maximum score
    max_score = 0
    chosen_test = 0
    for e in this_problem.U: #We compute the score for each test and update the maximum score test
        score = 0.0
        num_copies_positive = 0.0 #|{all copies of i, for i in H^+}|
        num_copies_negative = 0.0 #|{all copies of i, for i in H^-}|
        positive_probability = 0.0 #Pr(|H^+|)
        negative_probability = 0.0 #Pr(|H^-|)
        star_probability = 0.0 #Pr(|H^*|)
        compatible_Te_positive = this_problem.Te_positive[e] & this_problem.H
        compatible_Te_negative = this_problem.Te_negative[e] & this_problem.H
        positive_fraction = len(compatible_Te_positive)/(len(this_problem.H)-1) #This is the number of original compatible scnerios corresponding to + branch
        negative_fraction = len(compatible_Te_negative)/(len(this_problem.H)-1) #This is the number of original compatible scnerios corresponding to - branch
        for s in compatible_Te_positive:
            num_copies_positive = num_copies_positive + this_problem.num_copies[s]
        for s in compatible_Te_negative:
            num_copies_negative = num_copies_negative + this_problem.num_copies[s]
        for s in compatible_Te_positive:
            positive_probability = positive_probability + this_problem.p[s]
        for s in compatible_Te_negative:
            negative_probability = negative_probability + this_problem.p[s]
        for s in this_problem.H - (this_problem.Te_positive[e]|this_problem.Te_negative[e]):
            star_probability = star_probability + this_problem.p[s]
        if(num_copies_positive > num_copies_negative): #This if is for determining L_e side
            score = negative_probability * (positive_fraction + 1) + positive_probability * negative_fraction + 0.5 * star_probability * (1 - negative_fraction - positive_fraction)
        else:
            score = positive_probability * (negative_fraction + 1) + negative_probability * positive_fraction + 0.5 * star_probability * (1 - negative_fraction - positive_fraction)
        if(score > max_score):
            max_score = score
            chosen_test = e
    this_problem.U.discard(chosen_test) #We remove the selected element from U
    return chosen_test

def update_H(this_problem, e, realized_scenario): #Here we update the set of compatible scneraios (Scn) after showing test e
    feedback = this_problem.outcome[realized_scenario][e]
    if(feedback == 'u'): #If realized scenario has unknown outcome on e, we generate it here 
        feedback = random.getrandbits(1)
    incompatible_scenarios = set()
    for s in this_problem.H:
        if(this_problem.outcome[s][e] == 'u'): #Updating the probability of scenarios with unknown outcome for test e.
            this_problem.p[s] = this_problem.p[s]/2
            this_problem.num_copies[s] = this_problem.num_copies[s]/2
        elif(this_problem.outcome[s][e] != feedback): #Removing the incompatible scenarios
            incompatible_scenarios.add(s)
    this_problem.H = this_problem.H - incompatible_scenarios
    
def ODT(this_problem):
    cost = 0.0
    for realized_scenario in this_problem.H_all: #We run the algorithm for every scnerio as realized scenario
        this_problem.U = set()
        this_problem.H = set()         
        for e in this_problem.U_all:
            this_problem.U.add(e)
        for s in this_problem.H_all:
            this_problem.H.add(s)
            this_problem.p[s] = this_problem.pi[s]
            this_problem.num_copies[s] = this_problem.num_copies_all[s]
        while(len(this_problem.H) > 1): #Until there is more than one compatible scenario we continue 
            chosen_test = one_step_ODT(this_problem)
            update_H(this_problem, chosen_test, realized_scenario)
        cost = cost + (len(this_problem.U_all)-len(this_problem.U))*this_problem.pi[realized_scenario] #The cost of realized scneario is the number of discarded tests from U
    return cost

#Reading and initializing
file_name = input("Enter the data file name please: ")
rep = int(input("How many times do you want to run the algorithm? (Star outcomes are sampled everytime) ")) #Number of repetitions of algorithms (to have a reasonable sampling for star outcomes)

wb = open_workbook(file_name)
total_sheets = len(wb.sheet_names())

for i in range(total_sheets):
    wb_sheet = wb.sheet_by_index(i)
    this_problem = Problem(wb_sheet)
    this_problem.print_h_and_r()
    this_problem.compute_Te_positive()
    this_problem.compute_Te_negative()

    #Running the algorithm for each distribution for a number of repetitions and consider the average
    expected_cost = 0.0
    this_problem.init_probabilities(wb_sheet, alpha = 1)
    for i in range(rep):
        expected_cost = expected_cost + ODT(this_problem)  
    print("Algorithm cost for power law distribution with alpha = 0 (uniform distribution): ", expected_cost/rep)
    expected_cost = 0.0
    this_problem.init_probabilities(wb_sheet, alpha = 1.5)
    for i in range(rep):
        expected_cost = expected_cost + ODT(this_problem)    
    print("Algorithm cost for power law distribution with alpha = -0.5: ", expected_cost/rep)
    expected_cost = 0.0
    this_problem.init_probabilities(wb_sheet, alpha = 2)
    for i in range(rep):
        expected_cost = expected_cost + ODT(this_problem)    
    print("Algorithm cost for power law distribution with alpha = -1: ", expected_cost/rep)