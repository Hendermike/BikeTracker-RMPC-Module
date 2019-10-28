
public class RMPC {

	double Xprop;
	double Vprop;
	
	double Xant;
	double Vant;
	
	double Vlead;
	
	double desiredSpacing;
	
	int Npasos;
	double lb = -1.5;
	double ub = 1.0;
	double deltaLB = -3;
	double deltaUB = 1.5;
	
	//Restricciones
	double[][] Ades;
	double[] bdes;
	
	double[][] Aeq;
	double[] beq;
	
	//Modelo
	TSmodel model;
	
	//fuzzyInteral
	
	//Solo aguanta 5 pasos por ahora
	public RMPC(int Npasos, errorBuffer buffer) {

		//Ayuda memoria: pairStatePackage(double currentBicyclePositon, double currentBicycleSpeed, double previousBicyclePositon, double previousBicycleSpeed, double leaderSpeed)
		pairStatePackage paquete = new pairStatePackage(0, 5, 4, 5, 5);
		
		errorBuffer cloneBuffer = cloneBuffer(buffer);
		fuzzyInterval fuzzyInt = new fuzzyInterval();
		fuzzyInt.intervalosSim(cloneBuffer, Npasos);
		runGA(50, 0.60, 0.40, (int)(0.05*50), fuzzyInt, Npasos, paquete, 5000);
		
	}
	
	public void runGA(int populationSize, double mutationRate, double crossoverRate, int elitismCount, fuzzyInterval fuzzyInt, int Npasos, pairStatePackage paquete, int stopCondition) {
		
		GeneticAlgorithm ga = new GeneticAlgorithm(populationSize, mutationRate, crossoverRate, elitismCount, fuzzyInt);

		// Initialize population
		Population population = ga.initPopulation(2*Npasos);

		// Evaluate population
		ga.evalPopulation(population, paquete);

		// Keep track of current generation
		int generation = 1;

		/**
		 * Start the evolution loop
		 * 
		 * Every genetic algorithm problem has different criteria for finishing.
		 * In this case, we know what a perfect solution looks like (we don't
		 * always!), so our isTerminationConditionMet method is very
		 * straightforward: if there's a member of the population whose
		 * chromosome is all ones, we're done!
		 */
		int counter = 0;
		//int stopCondition = 1000;
		int toggle = 0;
		
		double tolerance = 0.000001;
		double[] fitArray = new double[] {0, 0};
		
		while (counter < stopCondition) {
			if (toggle == 0) {
				fitArray[toggle] = population.getFittest(0).getFitness();
				double delta = Math.abs(fitArray[0] - fitArray[1]);
				
				if (delta < tolerance) {counter += 1;}
				else {counter = 0;}
				
				toggle = 1;
				}
			else {
				fitArray[toggle] = population.getFittest(0).getFitness();
				double delta = Math.abs(fitArray[0] - fitArray[1]);
				
				if (delta < tolerance) {counter += 1;}
				else {counter = 0;}
				
				toggle = 0;
				}
			
			// Apply crossover
			population = ga.crossoverPopulation(population);

			// Apply mutation
			population = ga.mutatePopulation(population);

			// Evaluate population
			ga.evalPopulation(population, paquete);

			// Increment the current generation
			generation++;
		}

		/**
		 * We're out of the loop now, which means we have a perfect solution on
		 * our hands. Let's print it out to confirm that it is actually all
		 * ones, as promised.
		 */
		System.out.println("Found solution in " + generation + " generations");
		System.out.println("Best solution: " + population.getFittest(0).toString());
	}
	
	public errorBuffer cloneBuffer(errorBuffer buffer) {
		
		errorBuffer clone = new errorBuffer(buffer.getBufferLength());
		
		for (int i = 0; i < buffer.getBufferLength(); i++) {
			clone.setError(i, buffer.getError(i)); 
		}
		
		return clone;
	}
	
	public void defRestricciones(int Npasos, double lb, double ub, double deltaLB, double deltaUB) {
		double dt = 1;
		
		//Ades
		double[][] Alb = new double[2*Npasos][2*Npasos];
		double[][] Aub = new double[2*Npasos][2*Npasos];
		
		double[][] AdeltaLB = new double[Npasos - 1][Npasos];
		double[][] AdeltaUB = new double[Npasos - 1][Npasos];
		
		//Inicializando matrices
		for (int i = 0; i < 2*Npasos; i++) {
			for (int j = 0; j < 2*Npasos; j++) {
				if (i==j) {
					Alb[i][j] = -1;
					Aub[i][j] = 1;	
				} else {
					Alb[i][j] = 0;
					Aub[i][j] = 0;
				}
				if (i < Npasos - 1 && j < Npasos) {
					AdeltaLB[i][j] = 0;
					AdeltaUB[i][j] = 0;
				}
			}
		}
		
		//Filling delta matrix
		for (int i = 0; i < Npasos -1; i++) {
			AdeltaLB[i][i] = 1;
			AdeltaLB[i][i + 1] = -1;
			
			AdeltaUB[i][i] = -1;
			AdeltaUB[i][i + 1] = 1;
		}
		
		Ades = new double[8*Npasos - 4][6*Npasos];
		for (int i = 0; i < 8*Npasos - 4; i++) {
			for (int j = 0; j < 6*Npasos; j++) {
					if (i < 2*Npasos) {
						if(j < 4*Npasos) {Ades[i][j] = 0;}
						else {Ades[i][j] = Alb[i][j - 4*Npasos];}
					}
					else if (i < 4*Npasos && i >= 2*Npasos) {
						if(j < 4*Npasos) {Ades[i][j] = 0;}
						else {Ades[i][j] = Aub[i - 2*Npasos][j - 4*Npasos];}
					}
					else if (i < 5*Npasos - 1 && i >= 4*Npasos) {
						if(j < 4*Npasos) {Ades[i][j] = 0;}
						else if(j < 5*Npasos && j >= 4*Npasos) {Ades[i][j] = AdeltaLB[i - 4*Npasos][j - 4*Npasos];}
						else {Ades[i][j] = 0;}
					}
					else if (i < 6*Npasos - 2 && i >= 5*Npasos - 1) {
						if(j < 4*Npasos) {Ades[i][j] = 0;}
						else if(j < 5*Npasos && j >= 4*Npasos) {Ades[i][j] = 0;}
						else {Ades[i][j] = AdeltaLB[i - (5*Npasos - 1)][j - 5*Npasos];}
					}
					else if (i < 7*Npasos - 3 && i >= 6*Npasos - 2) {
						if(j < 4*Npasos) {Ades[i][j] = 0;}
						else if(j < 5*Npasos && j >= 4*Npasos) {Ades[i][j] = AdeltaUB[i - (6*Npasos - 2)][j - 4*Npasos];}
						else {Ades[i][j] = 0;}
					}
					else if (i < 8*Npasos - 4 && i >= 7*Npasos - 3) {
						if(j < 4*Npasos) {Ades[i][j] = 0;}
						else if(j < 5*Npasos && j >= 4*Npasos) {Ades[i][j] = 0;}
						else {Ades[i][j] = AdeltaUB[i - (7*Npasos - 3)][j - 5*Npasos];}
					}
			}
		}
		
		//bdes
		bdes = new double[8*Npasos - 4];
		for (int i = 0; i < 8*Npasos - 4; i++) {
			if (i < 2*Npasos) {bdes[i] = -lb;}
			else if (i < 4*Npasos && i >= 2*Npasos) {bdes[i] = ub;}
			else if (i < 6*Npasos - 2 && i >= 4*Npasos) {bdes[i] = -dt;}
			else {bdes[i] = dt;}
		}
		
	}
	
	
	//Este metodo deja de servir al implemntar dentro del funcional de costos las restricciones de igualdad.
	public void matIgualdades(int Npasos) {
		
		double dt = 2;
		
		double[][] bloque = new double[][] {{-1, -dt, 0, 0},
			{0, -1, 0, 0},
			{0, 0, -1, -dt},
			{0, 0, 0, -1}};
		
		double[][] A1 = new double[4*Npasos][4*Npasos];
		for( int i = 0; i < Npasos; i++) {
			A1[i][i] = 1;
		}																		
		int columnas;									
		for (int filas = 4; filas < 4*Npasos; filas = filas + 4) {
			columnas = filas - 4;
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					A1[filas + i][columnas + j] = bloque[i][j];
					//System.out.println("A1["+(filas + i)+"]["+(columnas + j)+"] :" + A1[filas + i][columnas + j]);
					}
				}
			}
		
		double[][] A2 = new double[4*Npasos][2*Npasos];
		for (int i = 0; i < 4*Npasos; i++) {
			for (int j = 0; j < 2*Npasos; j++) {
				if(i == 4*j + 1 || i == 4*(j - Npasos + 1) - 1) {A2[i][j] = -dt;}
				else {A2[i][j] = 0;}
			}
		}
		
		Aeq = new double[4*Npasos][6*Npasos];
		for (int i = 0; i < 4*Npasos; i++) {
			for (int j = 0; j < 6*Npasos; j++)
			{
				if (j < 4*Npasos) {Aeq[i][j] = A1[i][j];}
				else {Aeq[i][j] = A2[i][j - 4*Npasos];}
				System.out.println("Aeq[" + i + "][" + j + "] :" + Aeq[i][j]);
			}
		}
		
		//beq
		beq = new double[4*Npasos];
		for (int i = 0; i < 4*Npasos; i++) {
			if (i == 0) {beq[i] = Xprop + dt*Vprop;}
			else if (i == 1) {beq[i] = Vprop;}
			else if (i == 2) {beq[i] = Xant + dt*Vant;}
			else if (i == 3) {beq[i] = Vant;}
			else {beq[i] = 0;}
		}
											
	}
	
	public static void main(String [] args)	{

		int Npasos = 5;
		errorBuffer buffer = new errorBuffer(new double[] {0.01, -0.3, 0.2, 0.21, 0.01});
		RMPC rmpc = new RMPC(Npasos, buffer);
	
	}
	
}
