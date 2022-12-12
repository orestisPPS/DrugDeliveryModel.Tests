using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.Constitutive.Structural.Continuum;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class BrainInitialStructuralOnslyTest
    {
        const double Sc = 0.1;
        const double timeStep = 1; // in days
        const double totalTime = 10; // in days

        static double miNormal = 5; //KPa
        static double kappaNormal = 6.667; //Kpa
        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 201.74; //Kpa
        static int currentTimeStep = 0;
        static double lambda0 = 1;
        static Dictionary<double, double[]> Solution = new Dictionary<double, double[]>(); private static List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>();

		public BrainInitialStructuralOnslyTest()
		{
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        

        

        [Theory]
        [InlineData("../../../DataFiles/workingBrainComsolMesh.mphtxt")]
        public void StaticLinearTest(string fileName)
        {
            var equationModel = new ElasticBrainModelBuilder(fileName, Sc, miNormal, kappaNormal, miTumor, kappaTumor, timeStep, totalTime, lambda0);
            Dictionary<int, double> lambda = new Dictionary<int, double>(equationModel.Reader.ElementConnectivity.Count());
            foreach (var elem in equationModel.Reader.ElementConnectivity)
            {
                lambda.Add(elem.Key,  1d);
            }
            var model = new Model[] { EquationModels.ElasticBrainModelBuilder.CreateElasticModelFromComsolFile(equationModel.Reader,  lambda), };
            var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            var algebraicModel = new[] { solverFactory.BuildAlgebraicModel(model[0]), };
            var solver = new[] { solverFactory.BuildSolver(algebraicModel[0]), };
            var problem = new[] { new ProblemStructural(model[0], algebraicModel[0]), };
            var linearAnalyzer = new LinearAnalyzer(algebraicModel[0], solver[0], problem[0]);
            var analyzer = new StaticAnalyzer( algebraicModel[0], problem[0], linearAnalyzer);
            analyzer.Initialize();
            analyzer.Solve();
        }

        

        

    }
}
