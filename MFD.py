from ortools.linear_solver import pywraplp
from ortools.init import pywrapinit

def get_edge(x):      return (int(x.name().split('_')[1]), int(x.name().split('_')[2]))
def flatten(l):       return [item for sublist in l for item in sublist]
def head(x):          return int(x.name().split('_')[1])
def tail(x):          return int(x.name().split('_')[2])
def exists_edge(e,p): return len(list(filter(lambda x: tail(e)==x[0] and head(e)==x[1], p)))>=1

'''def flatten(l):
    if l==[]:
        return []
    else:
        if not isinstance(l[0], list):
            return [l[0]] + flatten(l[1:])
        else:
            return flatten(l[0]) + flatten(l[1:])'''


def PrintSolution(solution):
    print('Problem solved in ', solver.WallTime(), 'ms')
    print('Printing solution...')
    return


def Encode(k,A,F):

    edge_variables = []
    pi_variables   = []
    weights        = []

    for i in range(1,k+1): #DANGER the order of the variables in edge_varibles and pi_variables must be the same
        edge_variables += [ [solver.BoolVar('x_{}_{}_{}'.format(u,v,i))            for u in range(V) for v in range(V) if A[u][v]==1] ]
        pi_variables   += [ [solver.IntVar(0, infinity,'p_{}_{}_{}'.format(u,v,i)) for u in range(V) for v in range(V) if A[u][v]==1] ]
        weights        += [  solver.IntVar(1, w_max,'w_{}'.format(i)) ]

    #3a, 3b
    for i in range(k):
        solver.Add( sum( filter(lambda edge: head(edge)==source, edge_variables[i]) ) == 1 )
        solver.Add( sum( filter(lambda edge: tail(edge)==target, edge_variables[i]) ) == 1 )

    #3c
    for i in range(k):
        for v in range(1,V-1): #find all wedges u->v->w for a fixed v
            solver.Add( sum( filter(lambda edge: tail(edge)==v, edge_variables[i]) ) - sum( filter(lambda edge: head(edge)==v, edge_variables[i]) ) == 0 )

    #5a
    for e in range(E):
        solver.Add( sum( pi_variables[i][e] for i in range(len(pi_variables)) ) == F[head(pi_variables[0][e])][tail(pi_variables[0][e])] )

    #5b
    for i in range(k):
        for pi,x in tuple(zip(pi_variables[i],edge_variables[i])):
            solver.Add( pi <= w_max*x )

    #5c
    for i in range(k):
        for pi in pi_variables[i]:
            solver.Add( pi <= weights[i] )

    #5d
    for i in range(k):
        for pi,x in tuple(zip(pi_variables[i],edge_variables[i])):
            solver.Add( pi >= weights[i] - (1 - x)*w_max )

    #Add some objective function (we are just interested in deciding whether or not these constraints form a feasible region)
    solver.Minimize(sum(flatten(pi_variables)) + sum(flatten(edge_variables)))

    #Example of a subpath constraint: R=([[(1,3),(3,5)],[(0,1)]]), means that we have 2 paths, the first one is 1-3-5. the second path is just a single edge 0-1
    def EncodeSubpathConstraints(R):
        path_variables=[]
        for i in range(1,k+1):
            path_variables += [[ solver.BoolVar('r_{}_{}'.format(i,j)) for j in range(len(R)) ]]

        #[[r_1_0, r_1_1], [r_2_0, r_2_1], [r_3_0, r_3_1]] 
        #7a
        for i in range(k):
            for j in range(len(R)):
                solver.Add( sum( filter( lambda e: exists_edge(e,R[j]), edge_variables[i] )) >= len(R[j])*path_variables[i][j] )

        #7b
        for r in range(len(R)):
            solver.Add( sum( path_variables[i][r] for i in range(len(path_variables)) ) >= 1 )
        return

    def EncodeInexactFlow(F_low,F_high):
        #9a (to be exchanged with constraint 5a)
        for e in range(E):
            solver.Add( sum( pi_variables[i][e] for i in range(len(pi_variables)) ) <= F_high[head(pi_variables[0][e])][tail(pi_variables[0][e])] )
            solver.Add( sum( pi_variables[i][e] for i in range(len(pi_variables)) ) >= F_low[head(pi_variables[0][e])][tail(pi_variables[0][e])] )
        return

    #EncodeSubpathConstraints([[(1,3),(3,5)],[(0,1)]])

    print('Number of variables =', solver.NumVariables())
    print('Number of constraints =', solver.NumConstraints())

    return

'''
def BinarySearch():
    global solver
    global infinity
    L = 1
    R = min(E,f)

    while (L<=R):
        m = L + (R-L)//2

        solver = pywraplp.Solver.CreateSolver('SCIP')
        if not solver: exit(0)
        infinity = solver.infinity()

        Encode(k,A,F)
        
        status = solver.Solve()
        print("status:",status)

        if status == pywraplp.Solver.INFEASIBLE:
            R = m
        else:
            L = m+1

    print("Found decomposition of size",L)
    return status
'''


def LinearSearch(A,F):
    global solver
    global infinity
    k = len(list(filter(lambda f_sv: f_sv!=0, F[source]))) #k = |{v in V : f(s,v)>0}| trivial lower bound
    k=3
    while k <= min(E,f):

        solver = pywraplp.Solver.CreateSolver('SCIP') #creating a new solver in each iteration because I was having trouble with the clear() method... zz
        if not solver: exit(0)
        infinity = solver.infinity()

        Encode(k,A,F)
        
        status = solver.Solve()

        if status == pywraplp.Solver.INFEASIBLE:
            print("Infeasible with k="+str(k))
            k=k+1
        else:
            print("Found decomposition of size",k)
            break
    return status


def main():
    #input:  A DAG G=(V,E) and a flow function C, both given as |V|x|V| matrices A and F such that A[i][j]=1 iff (i,j) is an edge and F[i][j]=flow from node i to node j.
    #        In this program, we assume that all weights of the flow decomposition are non negative integers. Furthermore, only the matrix F is relevant, as edges that do
    #        not pass any flow are irrelevant for the minimum flow decomposition (they will never be part of any path of the decomposition).
    #output: A flow decomposition (P, w) such that |P| is minimized.
        
        # s a b c d t
    A = [[0,1,1,0,0,0],
         [0,0,1,1,0,0],
         [0,0,0,1,0,0],
         [0,0,0,0,1,1],
         [0,0,0,0,0,1],
         [0,0,0,0,0,0]]

    F = [[0,6,7,0,0,0],
         [0,0,2,4,0,0],
         [0,0,0,9,0,0],
         [0,0,0,0,6,7],
         [0,0,0,0,0,6],
         [0,0,0,0,0,0]]

    global V
    global E
    global f
    global w_max
    global source
    global target

    V = len(A)
    E = sum(flatten(A))
    f = sum(F[0])
    w_max = max(map(max,F))+1 #the largest flow value of any edge. since wi>=0 for any path i, wi is at most the highest value of flow that passes through any edge. (plus 1 because why not)
    source = 0
    target = V-1

    solution = LinearSearch(A,F)
    #solution = BinarySearch(A,F)

    PrintSolution(solution)


if __name__ == "__main__":
    main()
    print('Exiting now')