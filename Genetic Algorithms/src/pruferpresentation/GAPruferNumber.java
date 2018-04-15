package pruferpresentation;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import edgepresentation.MEF2Instances;
import filesinout.ReadFiles;
import operator.Crossover;
import operator.Evaluation;
import operator.InitializeChromosome;
import operator.Mutations;
import operator.Selection;
import random.MyRandom;
import structures.Individual;
import structures.Population;
public class GAPruferNumber{
	private static double m_rate = 0.05;
	private static final double p_rate = 0.5; // crossover rate
	private static  int generation = 500;
	private static int populationLength = 100; // default length of population
	private static Crossover crossover = new Crossover();
	private static Mutations mutation = new Mutations();
	private static Selection selection = new Selection();
	private static Evaluation evaluation = new Evaluation();
/*
 * type 6- small
2lin105-2x1","30kroB100-5x6","35kroB100-5x5","36eil101-6x6","42rat99-6x7","4berlin52-2x2","4eil51-2x2",
"4eil76-2x2","4pr76-2x2","6berlin52-2x3","6pr76-2x3","6st70-2x3","8berlin52-2x4","9eil101-3x3","9eil51-3x3","9eil76-3x3","9pr76-3x3
"2lin105-2x1",
type 4_large
4i200a","4i200h","4i200x1","4i200x2","4i200z","4i400a","4i400h","4i400x1","4i400x2","4i400z"
type3_large
6i300","6i350","6i400","6i450","6i500

 */


	public static void main(String[] args) {
		MFEAPruferNumber mfeaPruferNumber = new MFEAPruferNumber();
		String[] test = {"5i500-304","6i500","7i70-21","9pr439-3x3","2lin105-2x1","7i70-21","6i300","6i350","6i400","6i450"};
		for(int k= 0; k < test.length; k++){
		ReadFiles.clusterReadFiles("test/"+test[k]+".clt");
//	    ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/" +args[0]+ ".clt");
		int  num_vertex = ReadFiles.num_vertex;
//		 m_rate = 1.0/num_vertex;


//		Windows windows = new Windows();
//		windows.runWindow("The best ever");
		 String dirname = "Result/"+test[k];
	      File d1 = new File(dirname);
	      // Bay gio tao thu muc.
	      d1.mkdirs();

		 
//
	      System.out.println("-------------------------------------------------------------------------------");

		  for(int loop = 0; loop < 30; loop ++){
			MyRandom.setSeed(loop);
		    // write to files
		    PrintWriter pw = null;
			PrintWriter pw3 = null;
			
			
				try {
					pw = new PrintWriter(new FileWriter(dirname+"/"+test[k]+"-seed("+loop+").opt", true));
					pw3 = new PrintWriter(new FileWriter(dirname+"/"+test[k]+"_seed("+loop+").gen", true));
				} catch (IOException e) {
					e.printStackTrace();
				}
				pw3.println("Generations \t"+test[k]);
				
				 long start = System.currentTimeMillis();
			// Initialize population
		Population population = new Population();
		InitializeChromosome initializeChromosome = new InitializeChromosome();	
		population.population = initializeChromosome.initiaizePopulation(ReadFiles.clusters,populationLength, MyRandom.r);


		for (int i = 0; i <  generation; i++) {
			
			Population subPop = new Population(); // store new individuals are generated
			double[] popFitness = population.populationFitness1();


			for (int j = 0; j < populationLength / 2; j++) {
				Individual father = new  Individual();
				Individual mother = new  Individual();

				//use tournament to  select Individual
				father = population.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));
				mother = population.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));


				Individual offspring1 =  new Individual();
				Individual offspring2 =  new Individual();

				// new  random number  [0;1]
				double d =  MyRandom.r.nextDouble();
				// if random number >  crossover rate,  add  the individual to population without crossover
		     	if (d< p_rate) {
			//	 else do crossover, generate two offsprings, then  do mutation
		     		ArrayList<int[]> childs = new ArrayList<int[]>(); 
		     		childs = crossover.pruferNumberCrossover(father.getGene1(), mother.getGene1(),
			    			mfeaPruferNumber.calCulateGenLength(ReadFiles.clusters), MyRandom.r);
			    	offspring1.setGene( childs.get(0));
			    	offspring2.setGene( childs.get(1));

		     	}else{
						 offspring1.setGene(mutation.mutationPrufer1(ReadFiles.clusters,
									MyRandom.r,father.getGene1(), evaluation.getStartPoint(ReadFiles.clusters), m_rate));
						 offspring2.setGene(mutation.mutationPrufer1(ReadFiles.clusters,
									MyRandom.r,mother.getGene1(), evaluation.getStartPoint(ReadFiles.clusters), m_rate));
		     	
		     	}
				subPop.addIndiv(offspring1);
				subPop.addIndiv(offspring2);

			    }
			
			
			// get the index of the best Individual  of old Population
			double bestFiness = popFitness[0];
			int bestFinessIndex = 0;
			for(int t = 0; t < populationLength; t++ ){
			 if(popFitness[t] < bestFiness){
				 bestFiness = popFitness[t];
				  bestFinessIndex = t;
			 }

			}
			  pw3.println(i+"\t"+popFitness[bestFinessIndex]);



			// add the best Individual  of old Population to new Population
			 
			subPop.population.set(MyRandom.r.nextInt(populationLength),
					population.population.get(bestFinessIndex));
					population = subPop;
		}
		pw3.close();
		long end = System.currentTimeMillis();
		double[] popFitness = population.populationFitness1();
		double bestFiness = popFitness[0];
		int bestFinessIndex = 0;
		for(int t = 1; t < populationLength; t++ ){
			 if(popFitness[t] < bestFiness){
				 bestFiness = popFitness[t];
				  bestFinessIndex = t;
			 }
			}

	    // print to files
//		System.out.println(popFitness[bestFinessIndex]);
//		 pw.println(popFitness[bestFinessIndex]);	
		
		MEF2Instances mef = new MEF2Instances();
		 System.out.println("|"+loop+"\t|"+test[k]+"\t|"+popFitness[bestFinessIndex]+"\t|"+mef.getDateFromMillis(end - start)+"\t|");
		 pw.println("Name: "+test[k]);
		 pw.println("Seed: " + loop);
		 pw.println("Fitness: " +popFitness[bestFinessIndex]);
		 pw.println("Time: " + mef.getDateFromMillis(end - start));
//		 for(int m = 0; m < ReadFiles.num_vertex; m ++){
////			 pw.println();
////			 for( int n = 0; n < ReadFiles.num_vertex; n++){
////				 if( population.getIndividual(bestFinessIndex).getGene()[m][n] > 0){
////					 population.getIndividual(bestFinessIndex).getGene()[m][n] = 1.0f; }
////				 pw.print((int) population.getIndividual(bestFinessIndex).getGene()[m][n]+"\t");
////				
////			 }
//		 }
		 pw.close();

	}
		  }
}
}

