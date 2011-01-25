from pyfann import libfann as fann
class MyFANN():
    def __init__(self,xdat,ydat,idxs):
        if shape(xdat)[0] != shape(ydat)[0]:
            raise Exception('dimension mismatch b/w x, y')

        nt = len(xdat)
        
        ny = shape(ydat)[1]
        nx = shape(xdat)[1]

        num_input = nx;
        num_output = ny;
        num_layers = 3;
        num_neurons_hidden = 3;
        desired_error =  0.2;
        max_epochs =2000;
        epochs_between_reports = 1000;

        net = fann.neural_net()
        net.create_standard_array([num_layers, num_input, num_neurons_hidden, num_output]);

        net.set_activation_function_hidden( fann.SIGMOID_SYMMETRIC);
        net.set_activation_function_output( fann.SIGMOID_SYMMETRIC);
        
        t = fann.training_data()
        
        t.set_train_data(xdat,ydat)
        nt = net.train_on_data(t,max_epochs,epochs_between_reports,desired_error)
        out = net.save( "xor_float.net");

        print net.get_training_algorithm()
        raise Exception()

        fann.train_on_file( "xor.data", max_epochs, epochs_between_reports, desired_error);

        out = net.save( "xor_float.net");
        
        net.destroy();
        



from pyevolve import G1DList, GSimpleGA, Selectors, Scaling, DBAdapters
#from random import seed, randint, random
 
def eval_polynomial(x, *coefficients):
    result = 0
    for exponent, coeff in enumerate(coefficients):
        result += coeff*x**exponent
    return result 

def generate_fitness_function(sample_points):
    def fitness_function(chromosome):
        score = 0
        for point in sample_points:
            delta = abs(eval_polynomial(point[0], *chromosome) - point[1])
            score += delta
            score = -score
        return score
    return fitness_function
 
def run_pfit():
    # Generate a random polynomial, and generate sample points from it
    seed()
     
    source_polynomial = []
    for i in xrange(randint(1, 5)):
     source_polynomial.append(randint(-20,20))
     
    sample_points = []
    for i in xrange(20):
     n = randint(-100, 100)
     sample_points.append((n, eval_polynomial(n, *source_polynomial)))
     
    # Create the population
    genome = G1DList.G1DList(5)
    genome.evaluator.set(generate_fitness_function(sample_points))
    genome.setParams(rangemin=-50, rangemax=50)
     
    # Set up the engine
    ga = GSimpleGA.GSimpleGA(genome)
    ga.setPopulationSize(1000)
    ga.selector.set(Selectors.GRouletteWheel)
     
    # Change the scaling method
    pop = ga.getPopulation()
    pop.scaleMethod.set(Scaling.SigmaTruncScaling)
     
    # Start the algorithm, and print the results.
    ga.evolve(freq_stats=10)
    print(ga.bestIndividual())
    print("Source polynomial: " + repr(source_polynomial))
    print("Sample points: " + repr(sample_points))
