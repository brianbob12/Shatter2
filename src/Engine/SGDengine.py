from typing import Callable,Any
from tensorflow import Variable,Tensor,GradientTape
from tensorflow import math 
from tensorflow.keras.optimizers import Optimizer,SGD,Adam

class SGDengine:

  #solves f(x)=0
  @staticmethod
  def solve(f:Callable[[Tensor],Tensor],startingX:Any,learningRate:float,iterations:int)->Any:
    x=Variable(startingX)
    optimizer:Optimizer=Adam(learningRate)
    for i in range(iterations):
      with GradientTape() as g:
        g.watch([x])
        error=math.reduce_mean(math.abs(f(x)))
        grads=g.gradient(error,[x])
      optimizer.apply_gradients(zip(grads,[x])) 
    print(error)
    return x
    