package Cluster_Tree;
import java.util.ArrayList;
import java.util.Random;
import javax.swing.*;
public class GA{
	private static double m_rate = 0.2;
	private static final double p_rate = 0.8; // crossover rate
	private static  Random r = new Random();  
	private static  int generation = 10000;
	public static int populationLength = 100; // default length of population
	private static Crossover crossover = new Crossover();
	private static Mutations mutation = new Mutations();
	private static Selection selection = new Selection();
	
	private  static int num_vertex = ReadFiles.num_vertex;
	
	
	
	public static void main(String[] args) {
		
	
		JFrame gf = new JFrame();
		gf.setVisible(true);
		gf.setSize(800, 800);
		gf.setTitle("the best Individual");
		gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		gf.setVisible(true);
		// Initialize population
		Population population = new Population();
		for( int i = 0; i < populationLength; i++  ){
	        Individual individual = new Individual();
	        individual.inintilizeIndividual();
			population.addIndiv(individual);
		}
		
		System.out.print("---------------------------------------------> CLUSTER TREE <--------------------------------------------------\n ");
		System.out.println("Crossover rate : " + p_rate);
		System.out.println("Mutation  rate : " + m_rate);
		System.out.println("The number of generation : " +generation);
		
		
		
		for (int i = 0; i <  generation; i++) {
			Population subPop = new Population(); // store new individuals are generated
			double[] popFitness = population.populationFitness();
			
			
			for (int j = 0; j < Population.populationLength / 2; j++) {
				Individual father = new  Individual();
				Individual mother = new  Individual();
				
				//use tournament to  select Individual 
				father = population.getIndividual(selection.touranmentSelection(popFitness, 2, r));
				mother = population.getIndividual(selection.touranmentSelection(popFitness, 2, r));		
				
				Individual offspring1 =  new Individual();
				Individual offspring2 =  new Individual();
			    
				// new  random number  [0;1]
				double d = 0 + (1 - 0) * r.nextDouble();
				// if random number >  crossover rate,  add  the individual to population without crossover 
		     	if (d > p_rate) {
					subPop.addIndiv(father);
					subPop.addIndiv(mother);
			     } else {
				// else do crossover, generate two offsprings, then  do mutation
					offspring1.setGene(crossover.clusterCrossover(father.getGene(), mother.getGene(),
   		    		       num_vertex, ReadFiles.clusters));
					offspring2.setGene(crossover.clusterCrossover(father.getGene(), mother.getGene(),
	   		    		       num_vertex, ReadFiles.clusters));
					
					if( d < m_rate){
						 offspring1.setGene(mutation.mutationClusterTree(offspring1.getGene(),
								 num_vertex, ReadFiles.clusters));
						 offspring2.setGene(mutation.mutationClusterTree(offspring2.getGene(), 
								 num_vertex, ReadFiles.clusters));
					}
					
					
				subPop.addIndiv(offspring1);			 
				subPop.addIndiv(offspring2); 	
				
			}
			}
			// get the index of the best Individual  of old Population
			double bestFiness = popFitness[0];
			int bestFinessIndex = 0;
			
			for(int t = 1; t < populationLength; t++ ){
			 if(popFitness[t] < bestFiness){
				 bestFiness = popFitness[t];
				  bestFinessIndex = t;
			 }
			}
			// paint graph	
			Paint  p = new Paint();
			p.weightMatrix = population.getIndividual(bestFinessIndex).getGene();
			p.fitness = popFitness[bestFinessIndex];
			p.num_vertex = num_vertex;
		    gf.add(p);
		    gf.setVisible(true);
		    
			// add the best Individual  of old Population to new Population
			subPop.population.set(r.nextInt(Population.populationLength), 
					population.getIndividual(bestFinessIndex));	
					population = subPop;
		}
		// the best Individual at all
	   //	System.out.print("The best prufer number after " + generation + " generations: \n");
	}
}

