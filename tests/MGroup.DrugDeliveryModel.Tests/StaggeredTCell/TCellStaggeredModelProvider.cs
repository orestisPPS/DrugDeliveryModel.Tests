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
using MGroup.DrugDeliveryModel.Tests.PreliminaryModels;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using MGroup.MSolve.Numerics.Interpolation;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class TCellStaggeredModelProvider
    {

        public TCellModelProvider TCellModelProvider { get; }

        public Model[] model;



        //TODO put analysis time stepping where it belongs1 (pithanws sto Coupled7and9eqsSolution.cs h Coupled7and9eqsModel.cs)
        private GenericAnalyzerState[] analyzerStates, nlAnalyzerStates;
        private IParentAnalyzer[] parentAnalyzers;
        private IChildAnalyzer[] nlAnalyzers;
        private ISolver[] parentSolvers;


        public int CurrentTimeStep { get; set; }

        public GenericAnalyzerState[] AnalyzerStates => analyzerStates;
        public GenericAnalyzerState[] NLAnalyzerStates => nlAnalyzerStates;
        public IParentAnalyzer[] ParentAnalyzers => parentAnalyzers;
        public IChildAnalyzer[] NLAnalyzers => nlAnalyzers;
        public ISolver[] ParentSolvers => parentSolvers;
        public ComsolMeshReader Reader => reader;

        private ComsolMeshReader reader;

        private Dictionary<int, double[][]> SolidVelocityAtElementGaussPoints;

        private double timeStep;
        private double totalTime;

        private int incrementsPerStep;

        public TCellStaggeredModelProvider(TCellModelProvider tCellModelProvider, ComsolMeshReader comsolReader,
                                           Dictionary<int, double[][]> velocityAtElementGaussPoints,
                                           double timeStep, double totalTime, int incrementsPerStep)
        {
            TCellModelProvider = tCellModelProvider;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;


            analyzerStates = new GenericAnalyzerState[1];
            nlAnalyzerStates = new GenericAnalyzerState[1];
            parentAnalyzers = new IParentAnalyzer[1];
            nlAnalyzers = new IChildAnalyzer[1];
            parentSolvers = new ISolver[1];

            reader = comsolReader;

            
            this.SolidVelocityAtElementGaussPoints = velocityAtElementGaussPoints;
            
            this.timeStep = timeStep;
            this.totalTime  = totalTime;
            this.incrementsPerStep = incrementsPerStep;

            // intialize array ofm models1.
            model = new Model[1];
        }


        public void CreateModel(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            foreach (var element in reader.ElementConnectivity)
            {
                // SolidVelocityAtElementGaussPoints[elem.Key][0][0] = -5d;
                // SolidVelocityAtElementGaussPoints[elem.Key][0][1] = -5d;
                // SolidVelocityAtElementGaussPoints[elem.Key][0][2] = -5d;
                
                var elementNodes = element.Value.Item2;
                var elementGpVelocities = new double[1][];
                //elementGpVelocities[0] = GetVelocityVectorFromCoordinates(elementNodes);
                //elementGpVelocities[0] = GetVelocityVectorFromTime(elementNodes, CurrentTimeStep, timeStep);
                elementGpVelocities[0] = GetVelocityVectorFromTimeAndCoordinates(elementNodes, CurrentTimeStep, timeStep);
                SolidVelocityAtElementGaussPoints[element.Key] = elementGpVelocities;
                
            }
            
            model = new Model[1];
            
            model[0] = TCellModelProvider.GetModel();
            TCellModelProvider.AddBoundaryConditions(model[0]);
            (analyzers[0], solvers[0], nlAnalyzers[0]) = TCellModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep);

            for (int i = 0; i < analyzers.Length; i++)
            {
                analyzers[i].Initialize(true);
                if (analyzerStates[i] != null)
                {
                    analyzers[i].CurrentState = analyzerStates[i];
                }

                if (nlAnalyzerStates[i] != null)
                {
                    nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
                }
            }
        }

        public void CreateModelFirstTime(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            if (!(CurrentTimeStep == 0))
            {
                foreach (var element in reader.ElementConnectivity)
                {
                    var elementNodes = element.Value.Item2;
                    var elementGpVelocities = new double[1][];
                    //elementGpVelocities[0] = GetVelocityVectorFromCoordinates(elementNodes);
                    //elementGpVelocities[0] = GetVelocityVectorFromTime(elementNodes, CurrentTimeStep, timeStep);
                    elementGpVelocities[0] = GetVelocityVectorFromTimeAndCoordinates(elementNodes, CurrentTimeStep, timeStep);
                    SolidVelocityAtElementGaussPoints[element.Key] = elementGpVelocities;
                    
                    // SolidVelocityAtElementGaussPoints[elem.Key][0][0] = -5d;
                    // SolidVelocityAtElementGaussPoints[elem.Key][0][1] = -5d;
                    // SolidVelocityAtElementGaussPoints[elem.Key][0][2] = -5d;
                }
            }


            model = new Model[1];

            model[0] = TCellModelProvider.GetModel();
            TCellModelProvider.AddBoundaryConditions(model[0]);
            if (CurrentTimeStep == 0)
            {
                TCellModelProvider.AddInitialConditions(model[0]);
            }
            (analyzers[0], solvers[0], nlAnalyzers[0]) = TCellModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep);

            for (int i = 0; i < analyzers.Length; i++)
            {
                analyzers[i].Initialize(true);
                if (analyzerStates[i] != null)
                {
                    analyzers[i].CurrentState = analyzerStates[i];
                }

                if (nlAnalyzerStates[i] != null)
                {
                    nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
                }
            }
        }
        
        private double[] GetVelocityVectorFromCoordinates(Node[] elementNodes)
        {
            var interpolation = InterpolationTet4.UniqueInstance;
            var quadrature = TetrahedronQuadrature.Order1Point1;
            var shapeFunctionValues = interpolation.EvaluateFunctionsAt(quadrature.IntegrationPoints[0]);
            var gpCoordinates = new double[3]; //{ dphi_dksi, dphi_dheta, dphi_dzeta}
            for (var i1 = 0; i1 < shapeFunctionValues.Length; i1++)
            {
                gpCoordinates[0] += shapeFunctionValues[i1] * elementNodes[i1].X;
                gpCoordinates[1] += shapeFunctionValues[i1] * elementNodes[i1].Y;
                gpCoordinates[2] += shapeFunctionValues[i1] * elementNodes[i1].Z;
            }
            var spatiallyDistributedVelocityVector = new double[3];
            spatiallyDistributedVelocityVector[0] = 0d;
            spatiallyDistributedVelocityVector[1] = 0d;
            spatiallyDistributedVelocityVector[2] = - 50d * gpCoordinates[2];
            return spatiallyDistributedVelocityVector;
        }
        
        private double[] GetVelocityVectorFromTime(Node[] elementNodes, int currentTimeStep, double timeStep)
        {
            var interpolation = InterpolationTet4.UniqueInstance;
            var quadrature = TetrahedronQuadrature.Order1Point1;
            var shapeFunctionValues = interpolation.EvaluateFunctionsAt(quadrature.IntegrationPoints[0]);
            var gpCoordinates = new double[3]; //{ dphi_dksi, dphi_dheta, dphi_dzeta}
            for (var i1 = 0; i1 < shapeFunctionValues.Length; i1++)
            {
                gpCoordinates[0] += shapeFunctionValues[i1] * elementNodes[i1].X;
                gpCoordinates[1] += shapeFunctionValues[i1] * elementNodes[i1].Y;
                gpCoordinates[2] += shapeFunctionValues[i1] * elementNodes[i1].Z;
            }
            var timeDistributedVelocityVector = new double[3];
            timeDistributedVelocityVector[0] = 0d;
            timeDistributedVelocityVector[1] = 0d;
            timeDistributedVelocityVector[2] = - 500d * (currentTimeStep * timeStep);
            return timeDistributedVelocityVector;
        }
        
        private double[] GetVelocityVectorFromTimeAndCoordinates(Node[] elementNodes, int currentTimeStep, double timeStep)
        {
            var interpolation = InterpolationTet4.UniqueInstance;
            var quadrature = TetrahedronQuadrature.Order1Point1;
            var shapeFunctionValues = interpolation.EvaluateFunctionsAt(quadrature.IntegrationPoints[0]);
            var gpCoordinates = new double[3]; //{ dphi_dksi, dphi_dheta, dphi_dzeta}
            for (var i1 = 0; i1 < shapeFunctionValues.Length; i1++)
            {
                gpCoordinates[0] += shapeFunctionValues[i1] * elementNodes[i1].X;
                gpCoordinates[1] += shapeFunctionValues[i1] * elementNodes[i1].Y;
                gpCoordinates[2] += shapeFunctionValues[i1] * elementNodes[i1].Z;
            }
            var timeDistributedVelocityVector = new double[3];
            timeDistributedVelocityVector[0] = - 5000d * ((currentTimeStep * timeStep) * gpCoordinates[0]);
            timeDistributedVelocityVector[1] = - 5000d * ((currentTimeStep * timeStep) * gpCoordinates[1]);
            timeDistributedVelocityVector[2] = - 5000d * ((currentTimeStep * timeStep) * gpCoordinates[2]) ;
            return timeDistributedVelocityVector;
        }
    }
}


