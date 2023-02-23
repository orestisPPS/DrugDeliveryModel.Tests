using MGroup.MSolve.Discretization.Dofs;
using System;
using System.Collections.Generic;

namespace MGroup.DrugDeliveryModel.Tests.Commons;



public class CoupledProblemBoundaryConditionSet
{
        //Dictionary<TumorGrowthModelEquations, IEnumerable<List<>>
        public CoupledProblemBoundaryConditionSet(TumorGrowthTestCases testCase)
        {
                var a = new List<double>();
                
        }

        public CoupledProblemBoundaryConditionSet(List<(BoundaryConditionsUtility.BoundaryConditionCase, IDofType[], double[][], double[])> list)
        {
                
        }
}