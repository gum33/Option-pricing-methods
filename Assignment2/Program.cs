using System;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Differentiation;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;

namespace Assignment2
{
    //Class for the  Surface stochastic volatility inspired  parametrization
    class SSVIparametrisasion
    {
        //Parameters for SSVI that remain unchanged
        private double eta, gamma, rho, alpha, beta;


        //Parameters that change or need to be accessed
        public double S { get; private set; } //Price of underlying
        public double r { get; private set; } //risk free rate
        public double T { get; set; } //Maturity
        public double K { get; set; }//Strike price

        //Construct the surface
        public SSVIparametrisasion(double eta, double gamma, double rho, double alpha, double beta, double r, double S)
        {
            this.eta = eta;
            this.gamma = gamma;
            this.rho = rho;
            this.alpha = alpha;
            this.beta = beta;
            this.r = r;
            this.S = S;
        }


        //Calculate Log_Moneyness
        public double Log_moneyness(double K)
        {
            return Math.Log(K * Math.Exp(-r * (T)) / S);
        }

        //Implied volatility using wssvi
        public double SSVI_Impliedvolatility(double t, double k)

        {
            if (t < 0)
                throw new System.ArgumentException("Need t >= 0");

            return Math.Sqrt(Wssvi(t, k) / T); //Returns Sigma calculated from the ssvi surface
        }


        //Calculate value from the surface
        public double Wssvi(double t, double k)
        {
            double theta = Math.Pow(alpha, 2) * (Math.Exp(Math.Pow(beta, 2) * t) - 1);
            double phi = eta / (Math.Pow(theta, gamma) * Math.Pow(1 + theta, 1 - gamma));
            double wssvi = theta / 2 * (1 + rho * phi * k + Math.Sqrt(Math.Pow(phi * k + rho, 2) + 1 - Math.Pow(rho, 2)));
            return wssvi;
        }

        //Monte carlo method not using parallel loops to compare performance
        public double MonteCarloEuropeanCallSlow(double T, double K, int nrofpaths, int timesteps)
        {
            if (T <= 0 || K <= 0 || nrofpaths <= 0 || timesteps <= 0)
                throw new System.ArgumentException("Need T,K, number of paths and timesteps > 0");

            this.T = T;
            this.K = K;
            GeneratePath paths = new GeneratePath(this); //Generate S(t)
            double TotalCall = 0;
            double St = S;
            List<double> path;
            //Usi multi threading to run For loops in parallel making the simulation faster
            for(int i = 0; i< nrofpaths; ++i)
            {
                path = paths.generate_path(timesteps, T, true);
                St = path[path.Count - 1];
                TotalCall += Math.Max(St - K, 0);

            }

            return Math.Exp(-r * T) * TotalCall / nrofpaths;
        }
        //Monte carlo method not using parallel loops to compare performance
        public double MonteCarloEuropeanPutSlow(double T, double K, int nrofpaths, int timesteps)
        {
            if (T <= 0 || K <= 0 || nrofpaths <= 0 || timesteps <= 0)
                throw new System.ArgumentException("Need T,K, number of paths and timesteps > 0");
            //use put call parity to find put option price
            double TotalCall = MonteCarloEuropeanCallSlow(T, K, nrofpaths, timesteps) - S + K * Math.Exp(-r * (T));
            return TotalCall;

        }

        //Monte carlo methods implmented with parallel computing
        //European call option monte carlo
        public double MonteCarloEuropeanCall(double T, double K, int nrofpaths, int timesteps)
        {
            if(T <= 0 || K <= 0 || nrofpaths <= 0 || timesteps <=0 )
                throw new System.ArgumentException("Need T,K, number of paths and timesteps > 0");

            this.T = T;
            this.K = K;
            GeneratePath paths = new GeneratePath(this); //Generate S(t)
            double TotalCall = 0;
            double St = S;
            List<double> path;
            //Use multi threading to run For loops in parallel making the simulation faster
            Parallel.For(0, nrofpaths, i =>
            {
                path = paths.generate_path(timesteps, T, true);
                St = path[path.Count - 1];
                TotalCall += Math.Max(St - K, 0);
              
            });

            return Math.Exp(-r * T) * TotalCall / nrofpaths; //Take average and scale back in time
        }

        //European put option with monte carlo
        public double MonteCarloEuropeanPut(double T, double K, int nrofpaths, int timesteps)
        {
            if (T <= 0 || K <= 0 || nrofpaths <= 0 || timesteps <= 0)
                throw new System.ArgumentException("Need T,K, number of paths and timesteps > 0");

            double TotalCall = MonteCarloEuropeanCall(T, K, nrofpaths, timesteps) - S + K * Math.Exp(-r * (T));
            return TotalCall;

        }

        //Asian arithmetic call with monte carlo
        public double AsianArithmeticCallMonteCarlo(double T, double[] Tsteps, double K, int nrofpaths, int timesteps)
        {
            if (T<= 0 || K <= 0 || nrofpaths <= 0 || timesteps <= 0)
                throw new System.ArgumentException("Need T,K, number of paths and timesteps > 0");

            this.T = T; 
            this.K = K;
            double Calloption = 0, Stm = 0;
            GeneratePath paths = new GeneratePath(this);
            List<double> path;

            Parallel.For(0, nrofpaths, i =>
            {
                Stm = 0;
                path = paths.generate_path(timesteps, T);
                foreach (double Tm in Tsteps) //For each T1,...,Tm in T to average for asian options
                {
                    //Get the index of time step Tm
                        int index = Convert.ToInt32(Tm*timesteps);
                        Stm += path[index];
                    
                }
                
                double Savg = Stm / Tsteps.Length; //average over S(t1)+,...,+S(Tm)
                Calloption += Math.Max(Savg - K, 0); //Payoff function for asian arithmetic call function
             
            });
            return Math.Exp(-r * T) * Calloption /nrofpaths; //Take average and scale back in time
        }

        //Asian arithmetic put with monte carlo
        public double AsianArithmeticPutMonteCarlo(double T, double[] Tsteps, double K, int nrofpaths, int timesteps)
        {
            if (T <= 0 || K <= 0 || nrofpaths <= 0 || timesteps <= 0)
                throw new System.ArgumentException("Need T,K, number of paths and timesteps > 0");

            this.T = T;
            this.K = K;
            double Calloption = 0, Stm; //values for monte carlo method
            GeneratePath paths = new GeneratePath(this);
            List<double> path;
            Parallel.For(0, nrofpaths, i =>
            {
                Stm = 0.0;

                path = paths.generate_path(timesteps, T);
                foreach (double Tm in Tsteps) //For each T1,...,Tm in T to average for asian options
                {
                    int index = Convert.ToInt32(Tm * timesteps);
                    Stm += path[index];
                }
                double Savg = Stm / Tsteps.Length; //average over S(t1)+,...,+S(Tm)
                Calloption += Math.Max(Savg - K, 0); //Payoff function for asian arithmetic put function

            });
            return Math.Exp(-r * T) * Calloption /nrofpaths; //Take average and scale back in time
        }


        //Pricing look back option
        public double PricingLookbackOpt(double T, double K, int nrofpaths, int timesteps)
        {
            if (T <= 0 || K <= 0 || nrofpaths <= 0 || timesteps <= 0)
                throw new System.ArgumentException("Need T,K, number of paths and timesteps > 0");
            this.T = T;
            this.K = K;

            GeneratePath paths = new GeneratePath(this);
           
            double TotalCall = 0;
            double St = 0;
            List<double> path;
            Parallel.For(0, nrofpaths, i =>
            {
                path = paths.generate_path(timesteps, T, true);
                St = path[path.Count-1];

                TotalCall += St - path.Min(); //payoff function

            });
            return Math.Exp(-r * T)*TotalCall / nrofpaths; //Take average and scale back in time
        }
    }

    //Object that creates the path 
    class GeneratePath 
    {
        private NumericalDerivative numericalDerivative; //MathNets in built finite difference solver
        public SSVIparametrisasion ssvi;
        public double Sinf { get; set; } // Save the smallest St for use in Lookback
        public double T, K, S; // Parameters needed for path

        public GeneratePath(SSVIparametrisasion ssvi)
        {
            
            numericalDerivative = new MathNet.Numerics.Differentiation.NumericalDerivative(); //to approximate derivaties with the finite diff method.
            this.ssvi = ssvi; //The surface we are working with
            this.S = ssvi.S;
            this.T = ssvi.T;
        }

        //Returns the last timestep of the path
        //steps is is the value of steps per year, T is maturity time, lookback is used when there is no strike price, and t0 if we dont want to start at timestep 0
        public List<double> generate_path(double steps, double Tm, bool lookback = false, double t0 = 0)
        {
            if (T <= 0 || steps <= 0)
                throw new System.ArgumentException("Need T and timesteps > 0");

            List<double> path = new List<double>();
            path.Add(S);
            double dt = 1 / steps; //Size of each timejump so we get the number of steps each year as steps
            
            double Xt = Math.Log(ssvi.S);
            //if using lookback we find set K as S
            if (lookback)
                this.K = S;
            else
                this.K = ssvi.K;
            //Loop trough t1,t2,...,T where
            for (double t = t0; t<Tm; t+=dt)
            {
                if (lookback==true) //K value for lookback is min(S(t)) for t0,...,tcurrent
                {
                    this.K = path.Min();
                }
                double sigma = Dupire_sigma(this.T-t, Math.Exp(Xt));
            
                //get a sample from normal distribution
                double Zk = Normal.Sample(0, 1);
                Xt +=  (ssvi.r - 0.5 * sigma) * dt + Math.Sqrt(sigma * dt) * Zk; //update Xt

                path.Add(Math.Exp(Xt));                   
            }
            return path;
        }

        //Calculate Sigmadup using dupires formula from Lemma 2.2 (Gawlikowicz and Siska)
        private double Dupire_sigma(double t, double St)
        {

            //Log moneyness
            double k = Math.Log(this.K * Math.Exp(-ssvi.r * t) / St);
            
            double[] vars = { t, k };
            Func<double[], double> wtk = vars => ssvi.Wssvi(vars[0], vars[1]); //the wtk function that we need uses the SSVI surface
            double wprime_at_t = numericalDerivative.EvaluatePartialDerivative(wtk, vars, 0, 1);  //Evaluate partial derivative d/dt
            

            return wprime_at_t / G(vars, wtk);
        }

        //the g(t,k) in the Dupire formula
        private double G(double[] vars, Func<double[], double> wtk)
        {
            double wprime_at_k = numericalDerivative.EvaluatePartialDerivative(wtk, vars, 1, 1); //Evaluate first order derivatie at k
            double w2prime_at_k = numericalDerivative.EvaluatePartialDerivative(wtk, vars, 1, 2);//Evaluate second order deraivative at k

            //Equation (2.7) to find g(t,k)
              
            return Math.Pow(1 - vars[1] * wprime_at_k / (2.0 * wtk(vars)), 2) - Math.Pow(wprime_at_k, 2) / 4.0 * (1 / 4.0 + 1 / wtk(vars)) + w2prime_at_k / 2;
        }

    }
    
   
    class BlackScholesFormula
    {
        //Calculate european call option with the black scholes formula
        public static double CalculateCallOptionPrice(double S, double r, double sigma, double K, double T)
        {
            if (sigma <= 0 || T <= 0 || K <= 0 || S <= 0)
                    throw new System.ArgumentException("Need sigma > 0, T > 0, K > 0 and S > 0.");

            double d1 = (Math.Log(S / K) + (r + Math.Pow(sigma, 2) / 2) * T) / (sigma * Math.Sqrt(T));
            double N_of_d1 = MathNet.Numerics.Distributions.Normal.CDF(0, 1, d1); //Normal distributions value N(d1)
            double d2 = d1 - sigma * Math.Sqrt(T);
            double N_of_d2 = MathNet.Numerics.Distributions.Normal.CDF(0, 1, d2); //Normal distributions value N(d2)
            double callOption = S * N_of_d1 - K * Math.Exp(-r * T) * N_of_d2;
            return callOption;
        }
        
        //calculate european put option with the black scholes formula
        public static double CalculatePutOptionPrice(double S, double r, double sigma, double K, double T)
        {
            if (sigma <= 0 || T <= 0 || K <= 0 || S <= 0)
                throw new System.ArgumentException("Need sigma > 0, T > 0, K > 0 and S > 0.");
            double callOption = CalculateCallOptionPrice(S, r, sigma, K, T);
            double putOption = Math.Exp(-r * T) * K - S + callOption;
            return putOption;
        }
  
    }
    class Program
    {
        //Method that does everything needed in task 2.1 gives the values for the table
        static void doTask21(SSVIparametrisasion ssvi)
        {
            double sigma;

            double[] K = new double[] { 102.5313, 105.1271, 103.8212 }; //Ks for on the money 
            double[] T = new double[] { 1, 2, 1.5 };          //Ts for on the money
            Console.WriteLine("Task 2.1");
            Console.WriteLine("On the money prices");
            
            //Loop to get call option prices "on the money"
            for (int i = 0; i < K.Length; ++i)
            {
                ssvi.T = T[i];
                sigma = ssvi.SSVI_Impliedvolatility(T[i], 0);
                double call_price = (BlackScholesFormula.CalculateCallOptionPrice(ssvi.S, ssvi.r, sigma, K[i], T[i]));
                Console.WriteLine("For K = {0}, T = {1}, C(0, S) = {2}", K[i], T[i], Math.Round(call_price, 2));
            }
            double k;
            Console.WriteLine("\nOther prices");
            //Values for "other prices"
            K = new double[] { 80, 80, 90, 120, 90 };
            T = new double[] { 1, 2, 1, 1, 2 };
            //Loop that gets call option prices where k is not 0
            for (int i = 0; i < K.Length; ++i)
            {
                ssvi.T = T[i];
                k = ssvi.Log_moneyness(K[i]);
                
                sigma = ssvi.SSVI_Impliedvolatility(T[i], k);
                double price;
                if (i < 3) //call options
                {
                    price =  (BlackScholesFormula.CalculateCallOptionPrice(ssvi.S, ssvi.r, sigma, K[i], T[i]));
                    Console.WriteLine("For K = {0}, T = {1}, Call option price = {2}", K[i], T[i], Math.Round(price, 2));
                }
                else //put options
                {
                    price = (BlackScholesFormula.CalculatePutOptionPrice(ssvi.S, ssvi.r, sigma, K[i], T[i]));
                    Console.WriteLine("For K = {0}, T = {1}, Put option price = {2}", K[i], T[i], Math.Round(price, 2));

                } 
            }
        }

        //Runs monte carlo for the "other values" in Task 2.1 using the answer designed for Task 2.2
        static void doTask22(SSVIparametrisasion ssvi, int timesteps, int nrofpaths)
        {
            Console.WriteLine("\nTask 2.2 compared to task 2.1");
         
            double price = 0;

            double[] K = new double[] { 102.5313, 105.1271, 103.8212, 80, 80, 90, 120, 90 }; //Values to test
            double[] T = new double[] { 1, 2, 1.5 , 1, 2, 1, 1, 2 };       //Compare with Task 2.1
            Stopwatch sw = new Stopwatch(); //for measuring how long this simulation takes
            sw.Start();
            for (int i = 0; i < T.Length; ++i)
            {
                if (i < 6) //Call options
                {
                    price = ssvi.MonteCarloEuropeanCall(T[i], K[i], nrofpaths, timesteps);
                    Console.WriteLine("For K = {0}, T = {1}, Call option price = {2}", K[i], T[i], Math.Round(price, 2));
                }
                else  //Put options
                {
                    price = ssvi.MonteCarloEuropeanPut(T[i], K[i], nrofpaths, timesteps);
                    Console.WriteLine("For K = {0}, T = {1}, Put option price = {2}", K[i], T[i], Math.Round(price, 2));
                }
            }
            sw.Stop();

            Console.WriteLine("\nTime with parallel for loops elapsed = {0}", sw.Elapsed);
        }

        //Does monte carlo simulation to get asian arithmetic call options as described in task 2.4
        static void doTask24(SSVIparametrisasion ssvi, int timesteps, int nrofpaths)
        {
          
            Console.WriteLine("\nValues for Task 2.4");
            
            double[][] Tsteps = new double[][] 
            {
                new double[] { 0.1, 0.2, 0.3, 0.4, 0.5 },
                new double[] { 0.25, 0.50, 0.75, 1.00 },
                new double[] { 0.5,  1 }
            };
            double[] T = { 0.5, 1, 1.5 };
            double K = 100;
         
            for(int i = 0; i<T.Length; i++)
            {
                double MCprice = ssvi.AsianArithmeticCallMonteCarlo(T[i], Tsteps[i], K, nrofpaths, timesteps);
                Console.WriteLine("Asian arithmetic call option MC price: {0} for T: {1} ", Math.Round(MCprice,2), string.Join(", ",  Tsteps[i]));
            }
            
        }


        //Do monte carlo simulations for lookback options as described in task 2.5
        static void doTask25(SSVIparametrisasion ssvi, int timesteps, int nrofpaths)
        {
            
            Console.WriteLine("\nValues for Task 2.5");
            double[] T = new double[] { 0.5, 1, 1.5 };
            double K = 100;
        
            foreach (double t in T)
            {
                double MCprice = ssvi.PricingLookbackOpt(t, K, nrofpaths, timesteps);
                Console.WriteLine("Pricing lookback option MC price: {0} for T: {1} ", Math.Round(MCprice,2), t);
            }
        }

        //Compare performance of monte carlo wtih parallel for loops and regular loops
        static void doTaskExtra(SSVIparametrisasion ssvi, int timesteps, int nrofpaths)
        {
            Console.WriteLine("\nCompare speed of parallel Monte Carlo and regular Monte Carlo");
            double price; //option price
            double[] K = new double[] { 80, 80, 90, 120, 90 }; //Values to test
            double[] T = new double[] { 1, 2, 1, 1, 2 };       //Compare with Task 2.1
           
            for (int pths = 100; pths <= nrofpaths; pths *= 10)
            {
                Stopwatch sw = new Stopwatch(); //Object to time this method
                Console.WriteLine("\nRunning European option price monte carlo simulation with {0} paths using regular loops", pths);
                sw.Start();
                for (int i = 0; i < T.Length; ++i)
                {
                    if (i < 3) //Call options
                    {
                        price = ssvi.MonteCarloEuropeanCallSlow(T[i], K[i], pths, timesteps);
                    }
                    else  //Put options
                    {
                        price = ssvi.MonteCarloEuropeanPutSlow(T[i], K[i], pths, timesteps);
                    }
                  
                }
                sw.Stop();
                Console.WriteLine("\nTask 2.2 completed with number of paths {0}", pths);
                Console.WriteLine("Time with regular loops elapsed = {0}", sw.Elapsed);
                sw = new Stopwatch(); //Object to time this method
                Console.WriteLine("\nRunning European option price monte carlo simulation with {0} paths using parallel loops", pths);
                sw.Start();
                for (int i = 0; i < T.Length; ++i)
                {
                    if (i < 3) //Call options
                    {
                        price = ssvi.MonteCarloEuropeanCall(T[i], K[i], pths, timesteps);
                    }
                    else  //Put options
                    {
                        price = ssvi.MonteCarloEuropeanPut(T[i], K[i], pths, timesteps);
                    }
                }
                sw.Stop();
                Console.WriteLine("\nTask 2.2 completed with number of pths {0}", pths);
                Console.WriteLine("Time with parallel loops elapsed = {0}", sw.Elapsed);
            }
            

           
        }

        //Use main to call other methods
        static void Main(string[] args)
        {
            
            //Declare parameters
            double S = 100;
            double n = 0.2, rho = 0.3, gamma = 0.7, alpha = Math.Sqrt(0.1), beta = 1, r = 0.025; //constants given in assignment description
            SSVIparametrisasion ssvi = new SSVIparametrisasion(n, gamma, rho, alpha, beta, r, S); //create the surface
            int timesteps = 128, nrofpaths = 100000; //Change this to make code faster

            //Use static methods to complete task in a organized way
            //doTask21(ssvi);
            Console.WriteLine("\nMonte Carlo methods using {0} different paths with {1} time steps", nrofpaths, timesteps);

            doTask22(ssvi, timesteps, nrofpaths);
            doTask24(ssvi, timesteps, nrofpaths);
            doTask25(ssvi, timesteps, nrofpaths);

            doTaskExtra(ssvi, timesteps, nrofpaths);

            Console.WriteLine("\nProgram execution completed press any key to continue");
            Console.ReadKey();
        }
    }
}
