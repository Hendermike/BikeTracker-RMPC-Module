//TODO: incorporar restricciones a la variación de la aceleración y revisar código en general.

import java.util.Random;

/**
 * The GeneticAlgorithm class is our main abstraction for managing the
 * operations of the genetic algorithm. This class is meant to be
 * problem-specific, meaning that (for instance) the "calcFitness" method may
 * need to change from problem to problem.
 * 
 * This class concerns itself mostly with population-level operations, but also
 * problem-specific operations such as calculating fitness, testing for
 * termination criteria, and managing mutation and crossover operations (which
 * generally need to be problem-specific as well).
 * 
 * Generally, GeneticAlgorithm might be better suited as an abstract class or an
 * interface, rather than a concrete class as below. A GeneticAlgorithm
 * interface would require implementation of methods such as
 * "isTerminationConditionMet", "calcFitness", "mutatePopulation", etc, and a
 * concrete class would be defined to solve a particular problem domain. For
 * instance, the concrete class "TravelingSalesmanGeneticAlgorithm" would
 * implement the "GeneticAlgorithm" interface. This is not the approach we've
 * chosen, however, so that we can keep each chapter's examples as simple and
 * concrete as possible.
 * 
 * @author bkanber
 *
 */
public class GeneticAlgorithm {
	
	double lb = -1.5;
	double ub = 1.0;
	double deltaLB = -3;
	double deltaUB = 1.5;
	
	fuzzyInterval fuzzyInterval;
	
	private int populationSize;
	private double mutationRate;
	private double crossoverRate;
	private int elitismCount;

	public GeneticAlgorithm(int populationSize, double mutationRate, double crossoverRate, int elitismCount, fuzzyInterval fuzzyInterval) {
		this.populationSize = populationSize;
		this.mutationRate = mutationRate;
		this.crossoverRate = crossoverRate;
		this.elitismCount = elitismCount;
		this.fuzzyInterval = fuzzyInterval;
	}

	public Population initPopulation(int chromosomeLength) {
		// Initialize population
		Population population = new Population(this.populationSize, chromosomeLength);
		return population;
	}

	public double calcFitness(Individual individual, pairStatePackage paquete) {
		//Posiciones más recientes de la bicicleta 1:this 2:delantera
		double X_b1_k = paquete.getCurrentBicyclePositon();
		double V_b1_k = paquete.getCurrentBicycleSpeed();
		double X_b2_k = paquete.getPreviousBicyclePositon();
		double V_b2_k = paquete.getPreviousBicycleSpeed();
		double Vl = paquete.getLeaderSpeed();
		
		double ts = 0.5;
		
		//Ponderadores del funcional
		double Wx = 1, Wv = 1, Wu = 1; 
		
		//int Npasos = individual.getChromosomeLength();
		int Npasos = 5;
		double[] x = new double[4*Npasos];
		
		//Ayuda memoria
		// individual.getGene(0) = u_b1_k;
		// individual.getGene(1) = u_b1_k+1;
		// individual.getGene(2) = u_b1_k+2;
		// individual.getGene(3) = u_b1_k+3;
		// individual.getGene(4) = u_b1_k+4;
		// individual.getGene(5) = u_b2_k;
		// individual.getGene(6) = u_b2_k+1;
		// individual.getGene(7) = u_b2_k+2;
		// individual.getGene(8) = u_b2_k+3;
		// individual.getGene(9) = u_b2_k+4;
		
		x[0] = X_b1_k + ts*V_b1_k;
		x[1] = V_b1_k + ts*individual.getGene(0);
		x[2] = X_b2_k + ts*V_b2_k;
		x[3] = V_b2_k + ts*individual.getGene(5);
		
		//OJO
		for (int j = 1; j < Npasos; j++) {
			x[I(j)] = x[I(j) - 4] + ts*x[I(j) - 3];
			x[I(j) + 1] = x[I(j) - 3] + ts*individual.getGene(j);
			x[I(j) + 2] = x[I(j) - 2] + ts*x[I(j) - 1];
			x[I(j) + 3] = x[I(j) - 1] + ts*individual.getGene(Npasos + j);
		}
		
		double Jx = 0, Jv = 0, Ju = 0;
		double Ldes = 2, d = 2;
		double factX = Ldes + d;
		
		//Calcular score
		//Wx 17
		for (int i = 0; i < 4*Npasos - 3; i = i + 4) {
			Jx += (x[i] - x[i + 2] + factX)*(x[i] - x[i + 2] + factX);
		}
		//Wv 18
		for (int i = 1; i < 4*Npasos - 2; i = i + 4) {
			Jv += (x[i] - Vl)*(x[i] - Vl);
		}
		//Wu 4
		for (int i = 0; i < Npasos - 1; i++) {
			Ju += (individual.getGene(i + 1) - individual.getGene(i))*(individual.getGene(i + 1) - individual.getGene(i));
		}
		
		double fitness = (Wx*Jx + Wv*Jv + Wu*Ju);//(Npasos*(Wx + Wv + Wu)));

		// Store fitness
		individual.setFitness(fitness);

		return fitness;
	}
	
	public int I(int j) {
		return 4*j;
	}
	
	public void evalPopulation(Population population, pairStatePackage paquete) {
		double populationFitness = 0;

		// Loop over population evaluating individuals and suming population
		// fitness
		for (Individual individual : population.getIndividuals()) {
			populationFitness += calcFitness(individual, paquete);
		}

		population.setPopulationFitness(populationFitness);
	}

	public boolean isTerminationConditionMet(Population population) {
		for (Individual individual : population.getIndividuals()) {
			if (individual.getFitness() == 1) {
				return true;
			}
		}

		return false;
	}

	public Individual selectParent(Population population) {
		// Get individuals
		Individual individuals[] = population.getIndividuals();

		// Spin roulette wheel
		double populationFitness = population.getPopulationFitness();
		double rouletteWheelPosition = Math.random() * populationFitness;

		// Find parent
		double spinWheel = 0;
		for (Individual individual : individuals) {
			spinWheel += individual.getFitness();
			if (spinWheel >= rouletteWheelPosition) {
				return individual;
			}
		}
		return individuals[population.size() - 1];
	}

	//TODO: implementar crossover real (¿Es realmente necesario? GO...).
	public Population crossoverPopulation(Population population) {
		// Create new population
		Population newPopulation = new Population(population.size());

		// Loop over current population by fitness
		for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
			Individual parent1 = population.getFittest(populationIndex);

			// Apply crossover to this individual?
			if (populationIndex < this.elitismCount + (int)this.crossoverRate*population.size() && populationIndex >= this.elitismCount) {
				
				// Initialize offspring
				Individual offspring = new Individual(parent1.getChromosomeLength());
				
				// Find second parent
				Individual parent2 = selectParent(population);
				
				double alpha = Math.random();
				
				// Loop over genome
				for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {
					double newGene = (1 - alpha)*parent1.getGene(geneIndex) + alpha*parent2.getGene(geneIndex);
					offspring.setGene(geneIndex, newGene);
				}

				// Add offspring to new population
				newPopulation.setIndividual(populationIndex, offspring);
			} else {
				// Add individual to new population without applying crossover
				newPopulation.setIndividual(populationIndex, parent1);
			}
		}

		return newPopulation;
	}

	public Population mutatePopulation(Population population) {
		// Initialize new population
		Population newPopulation = new Population(this.populationSize);

		// Loop over current population by fitness: a partir de los no crossoveriados TODO
		for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
			Individual individual = population.getFittest(populationIndex);

			if (populationIndex >= this.elitismCount + (int)this.crossoverRate*population.size()) {
			// Loop over individual's genes
			for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
						if (Math.random() > 0.5) {
						double newGene = mutateGene(geneIndex, individual);
						// Mutate gene
						individual.setGene(geneIndex, newGene);
						}
				}
			}

			// Add individual to population
			newPopulation.setIndividual(populationIndex, individual);
		}

		// Return mutated population
		return newPopulation;
	}
	
	public double mutateGene(int geneIndex, Individual individual) {
		Random r = new Random();
		double gene = individual.getGene(geneIndex);
		
		double compressedUb = ub - fuzzyInterval.getSup(geneIndex%5);
		double compressedLb = lb + fuzzyInterval.getInf(geneIndex%5);
		
		double randomValue = r.nextGaussian()*(compressedUb - compressedLb)/12;
		
		gene = gene + randomValue;
		
		if (gene > compressedUb) {gene = compressedUb;}
		else if (gene < compressedLb) {gene = compressedLb;}
		
		//Segunda compresión
		return gene;
	}
	
	public void restrictionImplementation(Individual individual) {
		for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
			
			int Npasos = 5;
			
			double compressedUb = ub - fuzzyInterval.getSup(geneIndex%5);
			double compressedLb = lb + fuzzyInterval.getInf(geneIndex%5);
			
			double gene = individual.getGene(geneIndex);
			
			if (gene > compressedUb) {gene = compressedUb;}
			else if (gene < compressedLb) {gene = compressedLb;}
			
			if (geneIndex < individual.getChromosomeLength() - 1) {
				
				if (geneIndex != Npasos - 1) {
					double deltaU = individual.getGene((geneIndex + 1)) - gene;
					double adjuster;
					
					if (deltaU >= deltaUB) {
						adjuster = (deltaU - deltaUB);
					}
					else if (deltaU <= deltaLB) {
						adjuster = (deltaU - deltaLB);
					}
					else {adjuster = 0;}
					
					double adaptedGene1 = individual.getGene(geneIndex) + adjuster/2;
					double adaptedGene2 = individual.getGene(geneIndex + 1) - adjuster/2;
					individual.setGene(geneIndex, adaptedGene1);
					individual.setGene(geneIndex + 1, adaptedGene2);
				}				
				
			}
		}
	}
	
	public int getElitismCount() {
		return this.elitismCount;
	}

}
