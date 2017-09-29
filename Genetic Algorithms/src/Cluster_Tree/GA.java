package Cluster_Tree;
import java.util.ArrayList;
import java.util.Random;
public class GA {
	private static double m_rate = 0.1;
	private static final double p_rate = 0.8; // crossover rate
	static Random r = new Random(100);  
	private static  int generation = 1000;
	public static int populationLength = 100; // default length of population
	private static Crossover crossover = new Crossover();
	private static Mutations mutation = new Mutations();
	private static Selection selection = new Selection();
	private static ReadFiles readfiles = new ReadFiles();
	private Evaluation evaluation = new Evaluation();
	private InitializeChromosome initializeChromosome = new InitializeChromosome();
	private  static int num_vertex = ReadFiles.num_vertex;

	
	
	public static void main(String[] args) {
		//Population pop = new Population();
		//pop.initPopulation(); // initialize population
		
		
		// Initialize population
		Population population = new Population();
		for( int i = 0; i < populationLength; i++  ){
	        Individual individual = new Individual();
	        individual.inintilizeIndividual();
			population.addIndiv(individual);
		}
		
		
		System.out.print("---------------------------------------------> CLUSTER TREE <--------------------------------------------------\n ");
		//System.out.print("The best prufer number at first generation:\n ");
		//pop.getBestIndividual().printIndiv(); // the best Individual  at first generation
		System.out.println("Crossover rate : " + p_rate);
		System.out.println("Mutation  rate : " + m_rate);
		System.out.println("The number of generation : " +generation);
		Individual ind = new Individual(); 
		for (int i = 0; i <  generation; i++) {
			Population subPop = new Population(); // store new individuals are generated
			for (int j = 0; j < Population.populationLength / 2; j++) {
				Individual father = new  Individual();
				Individual mother = new  Individual();
				
				
				//father = .selectIndiv(2, r); //use tournament to  select Individual 
				//mother = pop.selectIndiv(2, r);
				father = population.getIndividual(0);
				mother = population.getIndividual(1);
				// test of mutation
				double[][] mothers = mutation.mutationClusterTree(mother.getGene(), num_vertex, ReadFiles.clusters);
			    Individual offsprings =  new Individual();
				// new  random number  [0;1]
				double d = 0 + (1 - 0) * r.nextDouble();
				// if random number >  crossover rate,  add  the individual to population without crossover 
				//if (d > p_rate) {
					//subPop.population.add(father);
					//subPop.population.add(mother);
				//} else {
				
					// else do crossover, generate two offsprings, then  do mutation
					double[][] newIndividual = crossover.clusterCrossover(father.getGene(), mother.getGene(), num_vertex, ReadFiles.clusters);
					//offsprings.get(0).mutation();
					//offsprings.get(1).mutation();
					subPop.population.add(offsprings);
					subPop.population.add(offsprings);
				//}
			}
			// add the best Individual  of old Population to new Population
			//ind = pop.getBestIndividual();
			//subPop.population.set(r.nextInt(Population.defaultPopLength), ind);
			//pop = subPop;
		}
		// the best Individual at all
		System.out.print("The best prufer number after " + generation + " generations: \n");
		//pop.getBestIndividual().printIndiv();
	}
}
/* 
 DIVIDE  tree to three group and  we have three group  inluding  vertices .
  */
		