
public class fuzzyInterval {

	double alpha;

	double Ypred;
	double interval;
	double sup;
	double inf;
	
	double[][] intervals;
	
	TSmodel model;
	
	public fuzzyInterval() {
		this.alpha = 0.67;
		this.model = new TSmodel();
	}
	
	public void newInterval(errorBuffer buffer) {
		this.Ypred = predictY(buffer);
		this.interval = interval(buffer);
		this.sup = Ypred + alpha*interval;
		this.inf = Ypred - alpha*interval;
	}
	
	public void intervalosSim(errorBuffer errorBuffer, int Npasos) {
		this.intervals = new double[Npasos][2];
		iIntUpdate(0, errorBuffer);
		for (int i = 0; i < Npasos; i++){
			errorBuffer.update(model.modelOutput(errorBuffer.getError(0), errorBuffer.getError(4)));
			iIntUpdate(i, errorBuffer);
			}
		}
	
	public void iIntUpdate(int i, errorBuffer errorBuffer) {
		newInterval(errorBuffer);
		intervals[i][0] = getSup();
		intervals[i][1] = getInf();
	}

	
	public void intervalosSimulacion(TSmodel model, int Npasos) {
		
	}
	
	public double predictY(errorBuffer buffer) {
		this.Ypred = model.modelOutput(buffer.getError(0), buffer.getError(4));
		return Ypred;
	}
		
	public double interval(errorBuffer buffer) {
		
		double out = incertezaUnaEntrada(buffer.getError(0), buffer.getError(4));
		return out;
		
	}
	
	public double incertezaUnaEntrada(double e1, double e5) {

		double[] betar = gradosActivacion(e1, e5);
		double[][] phiT = proyeccionEntrada(betar, new double[] {1, e1, e5});
		double [] Ir = calculoIr(phiT, model.getP(), model.getSigma());
		
		double out = 0;
		for (int i = 0; i < Ir.length; i++) {
			out = out + betar[i]*Ir[i];
		}
		return out;
		
	}
	
	public double [] gradosActivacion(double e1, double e5) {
		
		double[] activaciones = new double[5];
		
		double sum = 0;
		
		for (int i = 0; i<5; i++) {
			activaciones[i] = model.ruleActivation(i + 1, e1, e5);
			sum = activaciones[i] + sum;
		}
		
		for (int i=0; i<activaciones.length; i++) {
			activaciones[i] = activaciones[i] * (1/sum);
		}
		
		return activaciones;		
	}
	
	public double [][] proyeccionEntrada(double[] betar, double[] z) {
		
		int nIn = z.length;
		int nReglas = betar.length;
		
		double[][] proyDeEntradas = new double[nReglas][nIn];
		for (int i = 0; i < nReglas; i++) {
			for (int j = 0; j < nIn; j++) {
				proyDeEntradas[i][j] = betar[i]*z[j];
			}
		}
		
		return proyDeEntradas;
	};
	
	public double [] calculoIr(double[][] phiT, double[][][] P, double[] sigmas) {
		
		int nReglas = phiT.length;
		
		double[] Ir = new double[nReglas];
		
		double[] temp1 = new double[3];
		double temp2;
		double factor;
		
		for (int i = 0; i < nReglas; i++) {
			
			temp2 = 0;
			for (int j = 0; j < phiT[0].length; j++) {
				
				temp1[j] = 0;
				for (int k = 0; k < phiT[0].length; k++) {
					temp1[j] = phiT[i][k]*P[i][k][j] + temp1[j];	
				}
				
				temp2 = temp1[j]*phiT[i][j] + temp2;
			}
			
			factor = (1 + temp2)*(1 + temp2);
			Ir[i] = sigmas[i]*factor;
		}
		
		return Ir;
	}
	
	public double getSup(int i) {
		return intervals[i][0];
	}
	
	public double getInf(int i) {
		return intervals[i][1];
	}
	
	public double getSup() {
		return sup;
	}
	
	public double getInf() {
		return inf;
	}
	
	public double getPrediction() {
		return Ypred;
	}
	
	
	public static void main(String[] args) {
		
		System.out.println(9%5);
		
	}
	
}
