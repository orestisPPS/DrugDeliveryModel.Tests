using System.Collections.Generic;
using MGroup.Constitutive.Structural;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.DataStructures;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Solution;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Solvers.Direct;
using Xunit;
using modelBuilder = MGroup.DrugDeliveryModel.Tests.TemplateModel.Hexa8ContinuumNonLinearCantileverDynamicExample;
using System.Linq;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.Solvers.DofOrdering;
using MGroup.Solvers.DofOrdering.Reordering;
using MGroup.NumericalAnalyzers.Dynamic;
using System;
using MGroup.DrugDeliveryModel.Tests.EquationModels;

namespace MGroup.DrugDeliveryModel.Tests.TemplateModel
{
	public static class ValidateVelocity
	{   // ORIGIN: PlateRevisitedLessConstrainted.cs from Preliminary2_Copy repo
		// changes: updated to new msolve and for dynamic problem implementation

		//Edw epilegontai oi parametroi ts analushs
		public static int NR_steps = 1;
        private static double timestep=0.001; // those are overwritten for the periodic example
        private static double totalTime=0.100; // those are overwritten for the periodic example

		//Oi parametroi tou mondelou kai twn constrains  allazoun sto modelBuilder dld PlateRevisitedLessConstrainted

		const double Sc = 0.1;
		

		static double miNormal = 5; //KPa
		static double kappaNormal = 6.667; //Kpa
		static double miTumor = 22.44; //Kpa
		static double kappaTumor = 201.74; //Kpa
		static int currentTimeStep = 0;
		static double lambda0 = 1;
		static double density = 1;
		static StructuralDof loadedDof= StructuralDof.TranslationX;
		static double load_value =0.01;
		static int loadedNodeId;

		[Fact]
		private static void RunSuddenLoadTest()
		{
			var modelPath = "../../../DataFiles/workingTetMesh155.mphtxt";
			var equationModel = new MonophasicEquationModelDynamic(modelPath, Sc, miNormal, kappaNormal, miTumor, kappaTumor, timestep, totalTime, lambda0,
				density, loadedDof,load_value);

			Model model = equationModel.CreateElasticModelFromComsolFile();
			equationModel.AddBottomBCs(model, 0.1,0);
			equationModel.AddStaticNodalLoadsTopNodes(model, 0.1, 0);
			loadedNodeId = equationModel.loadedNodeId;

			double[] computedDisplacements = SolveModelDynamic(model);
			Assert.True(Utilities.AreDisplacementsSame(modelBuilder.GetExpectedDisplacementsSuddenLoad(), computedDisplacements, tolerance: 2E-3));
		}

		[Fact]
		private static void RunSuddenLoadWithInitialDisplacementsTest()
		{
			modelBuilder.monitoredDof = StructuralDof.TranslationX;
			Model model = modelBuilder.CreateModel();
			modelBuilder.AddStaticNodalLoads(model);
			modelBuilder.AddInitialConditionsDisplacements(model);
			double[] computedDisplacements = SolveModelDynamic(model);
			//TODO TO PRWTO VMA PREPEI NA PETIETAI KAI H SUGKRISI NA XEKINA ME TO DEFTERO
			Assert.True(Utilities.AreDisplacementsSame(modelBuilder.GetExpectedDisplacementsSuddenLoadAndInitialConditionsDisplacements(), computedDisplacements, tolerance: 1E-4));
		}




		[Fact]
		private static void RunTransientTestNoDelay()
		{
			modelBuilder.monitoredDof = StructuralDof.TranslationX;
			Model model = modelBuilder.CreateModel();
			modelBuilder.AddTransientLoadNoDelay(model);

			double[] computedDisplacements = SolveModelDynamic(model);
			Assert.True(Utilities.AreDisplacementsSame(modelBuilder.GetExpectedDisplacementsTransientLoadNoDelayADINA(), 
																computedDisplacements, tolerance: 1e-5));
		}

		[Fact]
		private static void RunTransientTestWithDelay()
		{
			modelBuilder.monitoredDof = StructuralDof.TranslationX;
			Model model = modelBuilder.CreateModel();
			modelBuilder.AddTransientLoadWithDelay(model);

			double[] computedDisplacements = SolveModelDynamic(model);
			Assert.True(Utilities.AreDisplacementsSame(modelBuilder.GetExpectedDisplacementsTransientLoadWithDelayADINA(),
																computedDisplacements, tolerance: 1e-5));
		}
			[Fact]
		private static void RunTransientTestPeriodic()
		{
			timestep = 0.0005; totalTime = 0.16;
			modelBuilder.monitoredDof = StructuralDof.TranslationX;
			Model model = modelBuilder.CreateModel();
			modelBuilder.AddPeriodicTransientLoad(model);

			double[] computedDisplacements = SolveModelDynamic(model);
			timestep = 0.0005;totalTime = 0.08;
			Assert.True(Utilities.AreDisplacementsSame(modelBuilder.GetExpectedDisplacementsPeriodicLoadADINA(),
																computedDisplacements, tolerance: 1e-5));
		}

		private static double[] SolveModelDynamic(Model model)
		{
			var solverFactory = new SuiteSparseSolver.Factory() { DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NodeMajorReordering()) };
			var algebraicModel = solverFactory.BuildAlgebraicModel(model);
			var solver = solverFactory.BuildSolver(algebraicModel);
			var problem = new ProblemStructural(model, algebraicModel);

			var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder( algebraicModel, solver, problem, numIncrements: NR_steps)
			{
				ResidualTolerance = 1E-10,
				MaxIterationsPerIncrement = 100,
				NumIterationsForMatrixRebuild = 1
			};
			var loadControlAnalyzer = loadControlAnalyzerBuilder.Build();


			//var staticAnalyzer = new StaticAnalyzer(model, algebraicModel, problem, loadControlAnalyzer);
			var dynamicAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder( algebraicModel, problem, loadControlAnalyzer,
				timeStep: timestep, totalTime: totalTime);
			dynamicAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration();
			NewmarkDynamicAnalyzer parentAnalyzer = dynamicAnalyzerBuilder.Build();
			parentAnalyzer.ResultStorage = new ImplicitIntegrationAnalyzerLog();



			//logs 
			var node_A = loadedNodeId;

			var emptyBCs = new List<INodalBoundaryCondition>();


			var loggerA = new TotalLoadsDisplacementsPerIncrementLog(model.NodesDictionary[node_A], loadedDof, emptyBCs, algebraicModel,
				$"hexaContinuumCantileverDynamicResults.txt");
			

			loadControlAnalyzer.IncrementalLog = loggerA;
			List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>();
			watchDofs.Add((model.NodesDictionary[node_A], loadedDof));
			loadControlAnalyzer.LogFactory = new LinearAnalyzerLogFactory(watchDofs, algebraicModel);
			
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			int totalNewmarkstepsNum = (int)Math.Truncate(totalTime / timestep);
			var totalDisplacementOverTime = new double[totalNewmarkstepsNum];
			for (int i1 = 0; i1 < totalNewmarkstepsNum; i1++)
            {
				var timeStepResultsLog = parentAnalyzer.ResultStorage.Logs[i1];
				totalDisplacementOverTime[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(node_A), loadedDof];
			}


			return totalDisplacementOverTime;



		}
	}

}
