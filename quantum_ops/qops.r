# Requirements:
# Qiskit is required for the functions which interact with it. Qiskit is called
# via the 'reticulate' library. If this library fails to find your python
# installation, you will need to manually specify it using use_python(). See
# below for an example.


library(R.utils)
library(R.oo)

library(reticulate)
#use_python("%LOCALAPPDATA%/Programs/Python/Python38/python.exe")


#--------------------------------------------------------------------------------
# Linear Algebra Functions
#--------------------------------------------------------------------------------
vec_norm <- function(v) sqrt(sum(v^2))
normalize <- function(v) v / vec_norm(v)

tensor <- function(...) {
  args <- list(...)
  k <- args[[1]]

  for(i in 2:length(args)) {
    k <- kronecker(k, args[[i]])
  }

  return(k)
}


#--------------------------------------------------------------------------------
# Quantum Primitives
#--------------------------------------------------------------------------------
bra <- function(...) matrix(as.complex(c(...)), nrow=1)
ket <- function(...) matrix(as.complex(c(...)), ncol=1)

# Converts a binary string/vector or decimal digit to a qubit representation. The
# binary string or vector has its most significant bit on the left.
as.qubit <- function(value, length=NULL) {
  dec_val <- 0
  len <- 0
  
  if (is.character(value)) { #binary string to qubit
    dec_val <- strtoi(value, base=2)
    len <- 2^nchar(value)
  }
  else if (is.vector(value) && is.numeric(value)) {
    if (length(value) == 1) { #decimal value to qubit
      dec_val <- value
      len <- 2^ceiling(max(1, log2(value+1)))
    }
    else { #binary integer vector to qubit
      dec_val <- Reduce(function(x,y) x*2+y, value) #convert binary to integer
      len <- 2^length(value)
    }
  }
  else {
    throw("as.qubit - Invalid argument: ", value)
  }
  
  if (!is.null(length)) {
    if (dec_val >= (2^length)) {
      throw("qubit of size ", length, " cannot hold value ", dec_val)
    }
    len <- length
  }
  
  state <- rep(0, len)
  state[dec_val+1] <- 1
  return(ket(state))
}


#--------------------------------------------------------------------------------
# Constants
#--------------------------------------------------------------------------------
q0 <- ket(1, 0)
q1 <- ket(0, 1)


#--------------------------------------------------------------------------------
# Gates
#--------------------------------------------------------------------------------
P0 <- ket(1, 0) %*% bra(1, 0) #projector, not a unitary gate
P1 <- ket(0, 1) %*% bra(0, 1) #projector, not a unitary gate
I <- function(size) diag(size)
H <- (1/sqrt(2)) * matrix(c(1, 1, 1, -1), nrow=2)
X <- matrix(c(0, 1, 1, 0), nrow=2)
Y <- matrix(c(0, -1i, 1i, 0), nrow=2)
Z <- matrix(c(1, 0, 0, -1), nrow=2)
S <- P <- matrix(c(1, 0, 0, 1i), nrow=2)
T <- matrix(c(1, 0, 0, exp((pi * 1i)/4)), nrow=2)
U1 <- function(phase) matrix(c(1, 0, 0, exp(phase * 1i)), nrow=2)
#CX <- CNOT <- matrix(c(1,0,0,0, 0,1,0,0, 0,0,0,1, 0,0,1,0), nrow=4)
#CX <- CNOT <- matrix(c(1,0,0,0, 0,0,0,1, 0,0,1,0, 0,1,0,0), nrow=4)
#CZ <- matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,-1), nrow=4)
#SWAP <- matrix(c(1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1), nrow=4)


# Less efficient, but more clear, method of constructing a controlled gate. See the actual
# implementation of CU.build for a more efficient method.
#   cond    - A projector (e.g. P0 or P1)
#   gate    - The gate to apply to the target qubit
#   control - The control qubit
#   target  - The target qubit
#   size    - The size of the system (e.g. 4 qubits)
# CU.build <- function(cond, gate, control, target, size) {
#   g <- 1
# 
#   for (i in rep(0, control-1)) { #i is unused
#     g <- g %x% I(2)
#   }
# 
#   g <- g %x% cond
#   
#   for (i in rep(0, target-control-1)) {
#     g <- g %x% I(2)
#   }
#   
#   g <- g %x% gate
#   
#   for (i in rep(0, size-target)) {
#     g <- g %x% I(2)
#   }
#   
#   return(g)
# }

# https://quantumcomputing.stackexchange.com/a/5192
# Builds a matrix that applies a gate to a target qubit in an n-qubit system
# if the control qubit meets a condition (e.g. is 1).
CU.build <- function(cond, gate, control, target, size) {
  stopifnot(all(dim(cond) == c(2, 2)))
  stopifnot(all(dim(gate) == c(2, 2)))
  stopifnot((control > 0) && (target > 1))
  stopifnot((control < size) && (target <= size))
  stopifnot((control < target) && (target > control))

  if (control > 1) {
    g <- I(2^(control-1)) %x% cond
  }
  else {
    g <- cond
  }
  
  if ((target-control) > 1) {
    g <- g %x% I(2^(target-control-1))
  }

  g <- g %x% gate

  if (target < size) {
    g <- g %x% I(2^(size-target))
  }

  return(g)
}

CU.build.rev <- function(cond, gate, control, target, size) {
  stopifnot(all(dim(gate) == c(2, 2)))
  stopifnot(all(dim(gate) == c(2, 2)))
  stopifnot((control > 1) && (target > 0))
  stopifnot((control <= size) && (target < size))
  stopifnot((control > target) && (target < control))

  if (control < size) {
    g <- I(2^(size-control)) %x% gate
  }
  else {
    g <- gate
  }
  
  if ((control-target) > 1) {
    g <- g %x% I(2^(control-target-1))
  }
  
  g <- g %x% cond
  
  if (target > 1) {
    g <- g %x% I(2^(target-1))
  }

  return(g)
}

CU.std <- function(gate, control, target, size) {
  g0 <- CU.build(P0, I(2), control, target, size)
  g1 <- CU.build(P1, gate, control, target, size)

  return(g0 + g1)
}

CU.rev <- function(gate, control, target, size) {
  g0 <- CU.build.rev(P0, I(2), control, target, size)
  g1 <- CU.build.rev(P1, gate, control, target, size)

  return(g0 + g1)
}

# This implementation treats the LSB as control, and so applies the
# reverse conditional gate in a typical scenario.
CX <- CNOT <- function(control, target, size) CU.rev(X, target, control, size)
CY <- function(control, target, size) CU.rev(Y, target, control, size)
CZ <- function(control, target, size) CU.rev(Z, target, control, size)


# Arguments to CU.std are modified to account for LSB being treated as the first bit.
SWAP <- function(bit1, bit2, size) {
  ctrl <- min(bit1, bit2)
  tgt <- max(bit1, bit2)
  CU.std(X, size-(tgt-1), size-(ctrl-1), size) %*% CU.rev(X, tgt, ctrl, size) %*% CU.std(X, size-(tgt-1), size-(ctrl-1), size)
}



#--------------------------------------------------------------------------------
# Probability Functions
#--------------------------------------------------------------------------------
probs <- function(q) abs(q)^2

# Probability of bit in position n being 1.
# Example: bit_prob(q, 4), with q being a 4-bit qubit, gives the probability of of q being 1xxx.
bit_prob <- function(q, bit) {
  stopifnot((bit > 0) && (bit <= log2(length(q))))

  prob <- 0

  # Number of partitions in the bitstring where bit n is consecutive 1's
  partitions <- length(q) / (2^bit)
  part_size <- 2^(bit-1) #size of each partition

  # Loop over each partition of consecutive 1's in position 'bit', and add
  # their probabilities to the total.
  for (part in seq(part_size, length(q)-part_size, 2*part_size)) {
    start <- part + 1
    end <- part + part_size
    prob <- prob + sum(probs(q[start:end]))
  }

  return(prob)
}

basis_probs <- function(q) {
  values <- as.list(probs(q))
  names(values) <- intToBin(0:(length(q)-1))
  return(values)
}

basis_amps <- function(q) {
  values <- as.list(abs(q))
  names(values) <- intToBin(0:(length(q)-1))
  return(values)
}

test <- function(q) runif(1) <= abs(q[2])^2


#--------------------------------------------------------------------------------
# Quantum Circuit Functions
#--------------------------------------------------------------------------------
qc.quantum_circuit <- function(num_qubits) {
  qc <- Object()
  qc$num_qubits <- num_qubits
  qc$instructions <- list()
  qc$statevector <- ket(1, rep(0, (2^num_qubits)-1))
  qc$compiled_gate <- I(2^num_qubits)
  
  # TODO: quantum gate functions as member functions? E.g.:
  # qc$h <- function(bit) qc.h(qc, bit)
  # or
  # qc$h <- function(bit) qc.apply_single_qubit_gate(qc, H, bit)

  return(qc)
}

qc.reset_state <- function(circuit) circuit$statevector <- ket(1, rep(0, (2^circuit$num_qubits)-1))

qc.compile <- function(circuit) {
  if (length(circuit$instructions) > 0)
    circuit$compiled_gate <- Reduce("%*%", rev(circuit$instructions))
}

qc.run <- function(circuit) circuit$statevector <- circuit$compiled_gate %*% circuit$statevector


# Applies a single-qubit gate to a specified qubit within a circuit.
qc.apply_single_qubit_gate <- function(circuit, gate, bit) {
  stopifnot((bit > 0) && (bit <= circuit$num_qubits))

  num_bits <- circuit$num_qubits
  final_gate <- I(2^max(num_bits - bit, 0))
  final_gate <- final_gate %x% gate
  final_gate <- final_gate %x% I(2^max(bit-1, 0))

  circuit$instructions[[length(circuit$instructions) + 1]] <- final_gate
}

# Applies a multi-qubit gate to the specified qubits within a circuit.
qc.apply_multi_qubit_gate <- function(circuit, gate, bit1, bit2) {
  stopifnot((bit1 > 0) && (bit1 <= circuit$num_qubits))
  stopifnot((bit2 > 0) && (bit2 <= circuit$num_qubits))

  final_gate <- gate(bit1, bit2, circuit$num_qubits)
  circuit$instructions[[length(circuit$instructions) + 1]] <- final_gate
}

qc.h <- function(circuit, bit, ...) {
  for (b in c(bit, ...)) {
    qc.apply_single_qubit_gate(circuit, H, b)
    circuit$program[[length(circuit$program) + 1]] <- list("H", b)
  }
}

qc.x <- function(circuit, bit, ...) {
  for (b in c(bit, ...)) {
    qc.apply_single_qubit_gate(circuit, X, b)
    circuit$program[[length(circuit$program) + 1]] <- list("X", b)
  }
}

qc.y <- function(circuit, bit, ...) {
  for (b in c(bit, ...)) {
    qc.apply_single_qubit_gate(circuit, Y, b)
    circuit$program[[length(circuit$program) + 1]] <- list("Y", b)
  }
}

qc.z <- function(circuit, bit, ...) {
  for (b in c(bit, ...)) {
    qc.apply_single_qubit_gate(circuit, Z, b)
    circuit$program[[length(circuit$program) + 1]] <- list("Z", b)
  }
}

qc.s <- function(circuit, bit, ...) {
  for (b in c(bit, ...)) {
    qc.apply_single_qubit_gate(circuit, S, b)
    circuit$program[[length(circuit$program) + 1]] <- list("S", b)
  }
}

qc.t <- function(circuit, bit, ...) {
  for (b in c(bit, ...)) {
    qc.apply_single_qubit_gate(circuit, T, b)
    circuit$program[[length(circuit$program) + 1]] <- list("T", b)
  }
}

qc.cx <- function(circuit, control, target) {
  qc.apply_multi_qubit_gate(circuit, CX, control, target)
  circuit$program[[length(circuit$program) + 1]] <- list("CX", control, target)
}

qc.cy <- function(circuit, control, target) {
  qc.apply_multi_qubit_gate(circuit, CY, control, target)
  circuit$program[[length(circuit$program) + 1]] <- list("CY", control, target)
}

qc.cz <- function(circuit, control, target) {
  qc.apply_multi_qubit_gate(circuit, CZ, control, target)
  circuit$program[[length(circuit$program) + 1]] <- list("CZ", control, target)
}

qc.swap <- function(circuit, bit1, bit2) {
  qc.apply_multi_qubit_gate(circuit, SWAP, bit1, bit2)
  circuit$program[[length(circuit$program) + 1]] <- list("SWAP", bit1, bit2)
}


#--------------------------------------------------------------------------------
# Quantum Program Functions
#--------------------------------------------------------------------------------

# A program starts with a 'qubits N' line that sets the number of qubits in the
# circuit. The following instructions are in the form 'gate i [j k ...]', where
# 'gate' is the gate name (e.g. 'H') and i, j, k, etc. are the qubits that the gate
# will act on. Gate names are case-insensitive. Comments and blank lines allowed.
# Comments start with '#'.
#
# The following is a 2-qubit program, which puts qubit 1 into superposition, then
# entangles it with qubit 2 via a conditional NOT:
# qubits 2
# h 1
# cx 1 2

qc.parse_opcode <- function(opcode) {
  switch(
    tolower(opcode),
    "h"    = qc.h,
    "x"    = qc.x,
    "not"  = qc.x,
    "y"    = qc.y,
    "z"    = qc.z,
    "s"    = qc.s,
    "t"    = qc.t,
    "cx"   = qc.cx,
    "cnot" = qc.cx,
    "cy"   = qc.cy,
    "cz"   = qc.cz,
    "swap" = qc.swap,
    throw("Invalid opcode: ", opcode)
  )
}

# Parses a program, as described above, and returns a quantum circuit implementing that program.
qc.parse_program <- function(text) {
  program <- list()

  lines <- unlist(strsplit(text, "\n"))
  lines <- Filter(function(line) !startsWith(line, "#") && (gsub("\\s", "", line) != ""), lines)

  qubits <- unlist(strsplit(lines[1], " "))
  if ((length(qubits) != 2) || (qubits[1] != "qubits")) {
    throw("Missing valid 'qubits N' declaration on first line")
  }
  program$num_qubits <- as.numeric(qubits[2])

  # Parse instructions
  for (line in lines[-1]) {
    split <- unlist(strsplit(line, " "))
    args <- Filter(function(x) x != "", split)

    if (length(args) == 0) {
      throw("Invalid line: ", line)
    }

    opcode <- args[1]
    operands <- as.numeric(args[-1])

    instr <- list(func = qc.parse_opcode(opcode), args = as.list(operands))
    program$instructions <- c(program$instructions, list(instr))
  }

  circuit <- qc.quantum_circuit(program$num_qubits)
  
  for (i in program$instructions) {
    args <- append(list(circuit), i$args)
    do.call(i$func, args)
  }
  
  qc.compile(circuit)

  return(circuit)
}


# Converts a quantum circuit to a qiskit quantum circuit object
qc.to_qiskit_ciruit <- function(circuit) {
  qk <- import("qiskit")
  
  #q <- qk$QuantumRegister(circuit$num_qubits)
  #c <- qk$ClassicalRegister(circuit$num_qubits)
  
  qc <- qk$QuantumCircuit(circuit$num_qubits)
  
  for (line in circuit$program) {
    func <- switch(
      tolower(line[[1]]),
      "h"    = qc$h,
      "x"    = qc$x,
      "not"  = qc$x,
      "y"    = qc$y,
      "z"    = qc$z,
      "s"    = qc$s,
      "t"    = qc$t,
      "cx"   = qc$cx,
      "cnot" = qc$cx,
      "cy"   = qc$cy,
      "cz"   = qc$cz,
      "swap" = qc$swap,
      throw("Invalid opcode: ", opcode)
    )
    
    # Subtract 1 from the register number (qiskit starts at 0), then call the function.
    do.call(func, as.list(as.integer(unlist(line[-1]) - 1)))
  }

  return(qc)
}

# Uses qiskit to generate the diagram
# https://qiskit.org/
qc.draw_to_file <- function(circuit, filename) {
  mpl <- import("matplotlib")
  qc <- qc.to_qiskit_ciruit(circuit)
  
  diagram <- qc$draw(output = "mpl")
  diagram$savefig(filename, bbox_inches = "tight")
  mpl$pyplot$close(diagram)
}

# Uses qiskit to generate the diagram
# https://qiskit.org/
qc.draw_to_text <- function(circuit) {
  qc <- qc.to_qiskit_ciruit(circuit)

  diagram <- qc$draw(output = "mpl")
  return(as.character(diagram))
}