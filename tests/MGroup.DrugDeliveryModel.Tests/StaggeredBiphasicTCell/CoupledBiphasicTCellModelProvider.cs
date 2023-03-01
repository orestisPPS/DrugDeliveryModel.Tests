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

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class CoupledBiphasicTCellModelProvider
    {
        public FluidPhaseModelProvider FluidPhaseModelProvider { get; }
        public SolidPhaseModelProvider SolidPhaseModelProvider { get; }
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

        private Dictionary<int, double> lambda;
        private Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints;
        private Dictionary<int, double[]> div_vs;

        private double timeStep;
        private double totalTime;

        private int incrementsPerStep;

        public CoupledBiphasicTCellModelProvider(FluidPhaseModelProvider fluidPhaseModelProvider,
                                                 SolidPhaseModelProvider solidPhaseProvider,
                                                    TCellModelProvider tCellModelProvider,
                                                 ComsolMeshReader comsolReader,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            Dictionary<int, double[]> div_vs, double timeStep, double totalTime, int incrementsPerStep)
        {
            SolidPhaseModelProvider = solidPhaseProvider;
            FluidPhaseModelProvider = fluidPhaseModelProvider;
            TCellModelProvider = tCellModelProvider;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;


            analyzerStates = new GenericAnalyzerState[3];
            nlAnalyzerStates = new GenericAnalyzerState[3];
            parentAnalyzers = new IParentAnalyzer[3];
            nlAnalyzers = new IChildAnalyzer[3];
            parentSolvers = new ISolver[3];

            reader = comsolReader;

            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            this.div_vs = div_vs;

            this.timeStep = timeStep;
            this.totalTime  = totalTime;
            this.incrementsPerStep = incrementsPerStep;

            // intialize array ofm models1.
            model = new Model[3];
        }


        public void CreateModel(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            //---------------------------------------
            // WARNING: do not initialize shared dictionarys because they have been passed by refernce in ewuationModel bilders.
            //---------------------------------------


            // update Shared quantities of Coupled model
            //foreach (var elem in reader.ElementConnectivity)
            //{ 
            //    lambda[elem.Key]= lambda0;
            //}
            foreach (var elem in reader.ElementConnectivity)
            {
                pressureTensorDivergenceAtElementGaussPoints[elem.Key] = ((ConvectionDiffusionElement3D)model[0].ElementsDictionary[elem.Key]).pressureTensorDivergenceAtGaussPoints;
            }
            foreach (var elem in reader.ElementConnectivity)
            {
                div_vs[elem.Key] = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocityDivergence;
            }
            
            model = new Model[3];
            
            //Create model for eq78 (fluid pressure)
            model[0] = FluidPhaseModelProvider.GetModel();
            FluidPhaseModelProvider.AddBoundaryConditions(model[0]);
            (analyzers[0], solvers[0], nlAnalyzers[0]) = FluidPhaseModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq9 (hyper-elastic material)
            model[1] = SolidPhaseModelProvider.GetModel();
            SolidPhaseModelProvider.AddBoundaryConditions(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = SolidPhaseModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);
            
            model[3] = TCellModelProvider.GetModel();
            TCellModelProvider.AddBoundaryConditions(model[3]);
            (analyzers[3], solvers[3], nlAnalyzers[3]) = TCellModelProvider.GetAppropriateSolverAnalyzerAndLog(model[3], timeStep, totalTime, CurrentTimeStep);

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
                foreach (var elem in reader.ElementConnectivity)
                {
                    pressureTensorDivergenceAtElementGaussPoints[elem.Key] = ((ConvectionDiffusionElement3D)model[0].ElementsDictionary[elem.Key]).pressureTensorDivergenceAtGaussPoints;
                }
                foreach (var elem in reader.ElementConnectivity)
                {
                    div_vs[elem.Key] = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocityDivergence;
                }
            }


            model = new Model[3];
            
            model[0] = FluidPhaseModelProvider.GetModel();
            FluidPhaseModelProvider.AddBoundaryConditions(model[0]);
            if(CurrentTimeStep==0)
            {
                //Eq78ModelProvider.AddEq78ModelInitialConditions(model[0]);
            }
            (analyzers[0], solvers[0], nlAnalyzers[0]) = FluidPhaseModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            model[1] = SolidPhaseModelProvider.GetModel();
            SolidPhaseModelProvider.AddBoundaryConditions(model[1]);
            
            model[3] = TCellModelProvider.GetModel();
            TCellModelProvider.AddBoundaryConditions(model[3]);
            if (CurrentTimeStep == 0)
            {
                TCellModelProvider.AddInitialConditions(model[3]);
            }
            
            (analyzers[1], solvers[1], nlAnalyzers[1]) = SolidPhaseModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

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

        public void SaveStateFromElements()
        {
            SolidPhaseModelProvider.SaveStateFromElements(model[1]);
        }
    }
}
