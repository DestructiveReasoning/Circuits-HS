module Circuits (
        parallel,   -- Adds two resistors in parallel
        para,       -- Adds a list of resistors in parallel
        roots,      -- Finds roots of quadratic equation
        todeg,      -- Converts radians to degrees
        torad,      -- Converts degrees to radians
        Phasor,     -- Complex Numbers
        CoordSyst,  -- Coordinate System of Complex Numbers
        addp,       -- Adds two Phasors
        mulp,       -- Multiplies two Phasors
        subp,       -- Subtracts two Phasors
        divp,       -- Divides two Phasors
        coordp,     -- Returns the coordinate system of a Phasor
        realp,      -- Returns the real component of a Phasor
        complexp,   -- Returns the complex component of a Phasor
        magp,       -- Returns the magnitude of a Phasor
        phasep,     -- Returns the phase of a Phasor
        phasepd,    -- Returns the phase of a Phasor (in degrees)
        conjp,      -- Returns the complex conjugate of a Phasor
        convertp,   -- Converts a Phasor to the other coordinate system
        parallelp,  -- Adds two Phasors in parallel
        parap,      -- Adds a list of Phasors in parallel
        resistor,   -- Returns a Phasor describing a resistor
        capacitor,  -- Returns a Phasor describing a capacitor
        inductor,    -- Returns a Phasor describing an inductor
        disp
) where

import Graphics.Gnuplot.Simple

data CoordSyst = Rectangular|Polar deriving (Eq,Show)
type Phasor a = (a,a, CoordSyst)

--instance Show (a,a,CoordSyst) where
--show (a,b,c)
--    | c == Rectangular = (Prelude.show a) ++ " + j" ++ (Prelude.show b)
--    | otherwise = (Prelude.show a) ++ "<" ++ (Prelude.show b)

todeg :: (Floating a) => a -> a
todeg x = x * 180 / pi

torad :: (Floating a) => a -> a
torad x = x * pi / 180

parallel :: Float -> Float -> Float
parallel r1 0 = 0
parallel 0 r2 = 0
parallel r1 r2 = 
    let sum = 1.0 / r1 + 1.0 / r2
    in 1.0 / sum

para :: [Float] -> Float
para (x:[]) = x
para (x:xs) = x `parallel` (para xs)

roots :: Float -> Float -> Float -> ((Float,Float),(Float,Float))
roots a b c = 
    let discriminant = b^2 - 4*a*c
    in 
        if discriminant >= 0 
        then (((-b + discriminant**0.5)/(2*a),0.0),((-b - discriminant**0.5)/(2*a),0.0))
        else 
            let discriminant' = -discriminant
            in
                ((-b/(2*a),(discriminant' ** 0.5)/(2*a)),(-b/(2*a),-(discriminant' ** 0.5)/(2*a)))

conjp :: (Floating a) => Phasor a -> Phasor a
conjp (a,b,c)
    | c == Rectangular = (a,-b,c)
    | otherwise = convertp $ conjp (convertp (a,b,c))

realp :: (Floating a) => Phasor a -> a
realp (a,b,s)
    | s == Rectangular = a
    | otherwise = realp $ convertp (a,b,s)

complexp :: (Floating a) => Phasor a -> a
complexp (a,b,s)
    | s == Rectangular = b
    | otherwise = complexp $ convertp (a,b,s)

magp :: (Floating a) => Phasor a -> a
magp (a,b,s)
    | s == Polar = a
    | otherwise = magp $ convertp (a,b,s)

phasep :: (Floating a) => Phasor a -> a
phasep (a,b,s)
    | s == Polar = b
    | otherwise = phasep $ convertp (a,b,s)

phasepd :: (Floating a) => Phasor a -> a
phasepd x = todeg $ phasep x

coordp :: (Floating a) => Phasor a -> CoordSyst
coordp (_,_,syst) = syst

convertp :: (Floating a) => Phasor a -> Phasor a
convertp (x,y,syst)
    | syst == Polar = (x * (cos y), x * (sin y),Rectangular)
    | syst == Rectangular = (sqrt (x**2 + y**2), atan(y/x),Polar)
    | otherwise = (0,0,Rectangular)

addp :: (Floating a) => Phasor a -> Phasor a -> Phasor a
addp u v
    | (((coordp u) == (coordp v)) && (coordp v == Rectangular)) = ((realp u) + (realp v), (complexp u) + (complexp v), Rectangular)
    | (coordp u == Rectangular) = addp u (convertp v)
    | (coordp v == Rectangular) = addp (convertp u) v
    | otherwise = convertp $ addp (convertp u) (convertp v)

subp :: (Floating a) => Phasor a -> Phasor a -> Phasor a
subp u v
    | (((coordp u) == (coordp v)) && (coordp v == Rectangular)) = ((realp u) - (realp v), (complexp u) - (complexp v), Rectangular)
    | (coordp u == Rectangular) = subp u (convertp v)
    | (coordp v == Rectangular) = subp (convertp u) v
    | otherwise = convertp $ subp (convertp u) (convertp v)

mulp :: (Floating a) => Phasor a -> Phasor a -> Phasor a
mulp u v
    | ((coordp u == coordp v) && coordp v == Polar) = ((magp u) * (magp v), (phasep u) + (phasep v), Polar)
    | (coordp u == Polar) = mulp u (convertp v)
    | (coordp v == Polar) = mulp (convertp u) v
    | otherwise = ((realp u) * (realp v) - (complexp u) * (complexp v), (realp u) * (complexp v) + (complexp u) * (realp v), Rectangular)

divp :: (Floating a) => Phasor a -> Phasor a -> Phasor a
divp u v
    | ((coordp u == coordp v) && coordp v == Polar) = ((magp u) / (magp v), (phasep u) - (phasep v), Polar)
    | (coordp u == Polar) = divp u (convertp v)
    | (coordp v == Polar) = divp (convertp u) v
    | otherwise = convertp $ divp (convertp u) (convertp v)

mulkp :: (Floating a, Ord a) => a -> Phasor a -> Phasor a
mulkp k (a,b,c)
    | c == Polar = if (k < 0) then (k*a, b + pi, c) else (k*a,b,c)
    | otherwise = convertp $ mulkp k (convertp (a,b,c))

parallelp :: (Floating a, Eq a) => Phasor a -> Phasor a -> Phasor a
parallelp a b = if (magp b /= 0.0) then divp (mulp a b) (addp a b) else (0.0,0.0,Rectangular)

parap :: (Floating a, Eq a) => [Phasor a] -> Phasor a
parap (x:[]) = x
parap (x:xs) = x `parallelp` (parap xs)

resistor :: (Floating a) => a -> Phasor a
resistor r = (r,0.0,Rectangular)

capacitor :: (Floating a) => a -> a -> Phasor a
capacitor c w = (0.0,-1/(w*c),Rectangular)

inductor :: (Floating a) => a -> a -> Phasor a
inductor l w = (0.0,l*w,Rectangular)

disp :: (Floating a, Show a) => Phasor a -> String
disp (a,b,c)
    | c == Rectangular = (Prelude.show a) ++ " + j" ++ (Prelude.show b)
    | otherwise = (Prelude.show a) ++ " < " ++ (Prelude.show (todeg b))
