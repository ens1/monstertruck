//! Finite outward intervals for algebraic boundary certification.

#[derive(Clone, Copy, Debug)]
pub(super) struct Interval {
    pub(super) low: f64,
    pub(super) high: f64,
}

impl Interval {
    pub(super) fn exact(value: f64) -> Option<Self> {
        value.is_finite().then_some(Self {
            low: value,
            high: value,
        })
    }

    fn new(low: f64, high: f64) -> Option<Self> {
        (low.is_finite() && high.is_finite() && low <= high).then_some(Self { low, high })
    }

    pub(super) fn point(self) -> Option<f64> { (self.low == self.high).then_some(self.low) }

    pub(super) fn neg(self) -> Self {
        Self {
            low: -self.high,
            high: -self.low,
        }
    }

    pub(super) fn add(self, other: Self) -> Option<Self> {
        if self.point() == Some(0.0) {
            return Some(other);
        }
        if other.point() == Some(0.0) {
            return Some(self);
        }
        if let (Some(a), Some(b)) = (self.point(), other.point()) {
            // Knuth's error-free sum retains exactly representable knot values.
            let sum = a + b;
            if !sum.is_finite() {
                return None;
            }
            let recovered = sum - a;
            let remainder = (a - (sum - recovered)) + (b - recovered);
            if remainder == 0.0 {
                return Self::exact(sum);
            }
        }
        Self::new(
            (self.low + other.low).next_down(),
            (self.high + other.high).next_up(),
        )
    }

    pub(super) fn sub(self, other: Self) -> Option<Self> { self.add(other.neg()) }

    pub(super) fn mul(self, other: Self) -> Option<Self> {
        if self.point() == Some(0.0) || other.point() == Some(0.0) {
            return Self::exact(0.0);
        }
        if self.point() == Some(1.0) {
            return Some(other);
        }
        if other.point() == Some(1.0) {
            return Some(self);
        }
        if self.point() == Some(-1.0) {
            return Some(other.neg());
        }
        if other.point() == Some(-1.0) {
            return Some(self.neg());
        }
        if let (Some(a), Some(b)) = (self.point(), other.point()) {
            let product = a * b;
            // Above 2^-960, the exact product's least possible nonzero residual
            // is representable. An underflowed FMA residual cannot claim exactness.
            if product.is_finite()
                && product.abs() >= 2.0_f64.powi(-960)
                && a.mul_add(b, -product) == 0.0
            {
                return Self::exact(product);
            }
        }
        let products = [
            self.low * other.low,
            self.low * other.high,
            self.high * other.low,
            self.high * other.high,
        ];
        if products.iter().any(|p| !p.is_finite()) {
            return None;
        }
        Self::new(
            products
                .into_iter()
                .fold(f64::INFINITY, f64::min)
                .next_down(),
            products
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        )
    }

    pub(super) fn div(self, other: Self) -> Option<Self> {
        if other.low <= 0.0 && other.high >= 0.0 {
            return None;
        }
        if self.point() == Some(0.0) {
            return Self::exact(0.0);
        }
        if other.point() == Some(1.0) {
            return Some(self);
        }
        if other.point() == Some(-1.0) {
            return Some(self.neg());
        }
        if let (Some(a), Some(b)) = (self.point(), other.point()) {
            let quotient = a / b;
            let product = quotient * b;
            if quotient.is_finite()
                && product.is_finite()
                && product.abs() >= 2.0_f64.powi(-960)
                && (-quotient).mul_add(b, a) == 0.0
            {
                return Self::exact(quotient);
            }
        }
        let quotients = [
            self.low / other.low,
            self.low / other.high,
            self.high / other.low,
            self.high / other.high,
        ];
        if quotients.iter().any(|p| !p.is_finite()) {
            return None;
        }
        Self::new(
            quotients
                .into_iter()
                .fold(f64::INFINITY, f64::min)
                .next_down(),
            quotients
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        )
    }

    pub(super) fn sqrt(self) -> Option<Self> {
        if self.high < 0.0 {
            return None;
        }
        Self::new(
            self.low.max(0.0).sqrt().next_down().max(0.0),
            self.high.sqrt().next_up(),
        )
    }

    pub(super) fn magnitude(self) -> f64 { self.low.abs().max(self.high.abs()) }

    pub(super) fn midpoint(self) -> Option<f64> {
        if let Some(point) = self.point() {
            return Some(point);
        }
        let result = self.low * 0.5 + self.high * 0.5;
        (result.is_finite() && result >= self.low && result <= self.high).then_some(result)
    }
}
